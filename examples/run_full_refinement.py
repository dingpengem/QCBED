"""Edit these settings; use --quick to validate the full real-data pipeline."""

from dataclasses import replace
from datetime import datetime
from pathlib import Path
import argparse
import json

from qcbed import ForwardExecutionConfig
from qcbed.refinement_subset import select_refinement_categories
from qcbed.refinement_bounds import RefinementBounds, RefinementBoundRules, resolve_refinement_bounds
from qcbed.full_runner import (
    FullRefinementConfig,
    run_full_refinement,
    load_refinement_checkpoint,
)
from qcbed.refinement_driver import (
    parse_legacy_refinement,
    LegacyRefinementConfig,
)
from qcbed.refinement import build_refinement_parameter_layout
from qcbed.experimental import parse_experimental_intensities
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records


# ============================================================
# PATHS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = ROOT / "reference_data" / "FOR007.DAT"
EXPERIMENTAL_FILE = ROOT / "reference_data" / "in.dat"


# ============================================================
# REFINEMENT SETTINGS
# ============================================================

# Current reference data contain zero-variance points only where
# weight == 0, so use the explicit non-legacy masked objective.
OBJECTIVE_MODE = "masked"

# Invalid numerical proposals are rejected with +inf rather than
# terminating a production refinement run.
INVALID_TRIAL_POLICY = "reject"

# ------------------------------------------------------------
# Active refinement categories
# ------------------------------------------------------------
#
# First controlled production run:
#
#   1 × depth
#   6 × scan_coordinate
#
# Everything else remains fixed at its starting value:
#
#   absorption
#   debye_waller
#   reflection_multiplier
#
# Total active parameters = 7.
#
REFINEMENT_CATEGORIES = (
    "depth",
    "scan_coordinate",
)

# User-selected production constraints, NOT source-derived legacy limits.
BOUNDS_POLICY = "reject"
DEPTH_MIN_ANGSTROM = 600.0
DEPTH_MAX_ANGSTROM = 1200.0
SCAN_RELATIVE_BOUND = 0.20
# Dimensionless undilated fractional reciprocal-cell coordinate half-width.
# Actual half-width = max(abs(start)*SCAN_RELATIVE_BOUND, this explicit floor).
SCAN_ABSOLUTE_FALLBACK = 0.05

# Concise user-selected rules; all active targets expand automatically.
SPECIFIC_BOUNDS = {("depth", (0,)): (DEPTH_MIN_ANGSTROM, DEPTH_MAX_ANGSTROM)}
CATEGORY_BOUNDS = {
    "absorption": (0.5, 1.5),
    "debye_waller": (0.5, 1.5),
    "reflection_multiplier": (0.5, 1.5),
}
CATEGORY_RELATIVE_BOUNDS = {"scan_coordinate": SCAN_RELATIVE_BOUND}
CATEGORY_RELATIVE_FALLBACK = {"scan_coordinate": SCAN_ABSOLUTE_FALLBACK}


# ============================================================
# FORWARD EXECUTION
# ============================================================

# Validated batched CPU backend.
EXECUTION_BACKEND = "batched"

# Real reference benchmark showed batch 64 performs very well.
BATCH_SIZE = 64


# ============================================================
# OPTIMIZER SETTINGS
# ============================================================

# Source-backed hybrid:
#
#   ANNEAL -> AMOEBA
#
OPTIMIZER_MODE = "ANSI"

LOOPS = 1
SEED = 2026


# ============================================================
# LEGACY OPTIMIZER LIMITS
# ============================================================

# ANNEAL source MAXITS.
MAXFSIMP = 20

# PARAMIN source cycle limit.
MAXFPARA = 3

# AMOEBA evaluation limit.
#
# ANNEAL and AMOEBA share the same legacy MAXFSIMP field in the
# current driver architecture, so for ANSI these values must agree.
AMOEBA_MAX_EVALUATIONS = 20

# SIMPARA evaluation limit.
SIMPARA_MAX_EVALUATIONS = 20


# ============================================================
# PROGRESS / CHECKPOINT SETTINGS
# ============================================================

# Save a CBED snapshot every N objective evaluations.
SNAPSHOT_EVERY = 25

# Save a checkpoint every N objective evaluations.
CHECKPOINT_EVERY = 25


# ============================================================
# OUTER SAFETY LIMITS
# ============================================================

# First controlled production run.
MAX_TOTAL_EVALUATIONS = 150

# 30 minutes.
MAX_WALL_TIME = 1800.0


# ============================================================
# CONFIGURATION
# ============================================================

def make_config(original, *, quick=False, refinement_bounds=RefinementBounds(), bound_rules=None):
    cfg = replace(
        original,
        method=OPTIMIZER_MODE,
        loops=LOOPS,
        seed=SEED,
        maxfpara=MAXFPARA,
        execution=ForwardExecutionConfig(
            EXECUTION_BACKEND,
            BATCH_SIZE,
        ),
    )

    # These stages share the same legacy MAXFSIMP field.
    limits = {
        "ANNEAL": MAXFSIMP,
        "AMOEBA": AMOEBA_MAX_EVALUATIONS,
        "SIMPARA": SIMPARA_MAX_EVALUATIONS,
    }

    selected = {
        limits[stage]
        for stage in cfg.stages
        if stage in limits
    }

    if len(selected) > 1:
        raise ValueError(
            "selected stages share source MAXFSIMP; "
            "their configured limits must agree"
        )

    cfg = replace(
        cfg,
        maxfsimp=next(iter(selected), MAXFSIMP),
    )

    if quick:
        # Short real-data pipeline validation.
        #
        # This does not claim scientific convergence.
        cfg = replace(
            original,
            method="ANNE",
            loops=1,
            seed=SEED,
            maxfsimp=0,
            maxtemp=0.001,
            mintemp=0.0001,
            coolrate=0.5,
            execution=ForwardExecutionConfig(
                "batched",
                64,
            ),
        )

    return FullRefinementConfig(
        cfg,
        invalid_trial_policy=INVALID_TRIAL_POLICY,
        bounds_policy=BOUNDS_POLICY,
        refinement_bounds=refinement_bounds,
        bound_rules=bound_rules,
        refinement_categories=REFINEMENT_CATEGORIES,
        objective_mode="masked" if quick else OBJECTIVE_MODE,
        snapshot_every=None if quick else SNAPSHOT_EVERY,
        checkpoint_every=1 if quick else CHECKPOINT_EVERY,
        max_total_evaluations=20 if quick else MAX_TOTAL_EVALUATIONS,
        max_wall_time=600.0 if quick else MAX_WALL_TIME,
    )


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = parser.add_mutually_exclusive_group()

    mode.add_argument(
        "--quick",
        action="store_true",
        help="run a very short real-data pipeline validation",
    )

    mode.add_argument(
        "--resume",
        type=Path,
        help="resume from an existing refinement run directory",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="parse and print configuration without evaluating the objective",
    )

    parser.add_argument(
        "--max-total-evaluations",
        type=int,
        help=(
            "override cumulative objective-evaluation budget, "
            "including on resume"
        ),
    )

    parser.add_argument(
        "--max-wall-time",
        type=float,
        help=(
            "override cumulative wall-time budget in seconds, "
            "including on resume"
        ),
    )

    args = parser.parse_args()

    source = INPUT_FILE
    exp_path = EXPERIMENTAL_FILE

    # ========================================================
    # RESUME EXISTING RUN
    # ========================================================

    if args.resume:
        checkpoint = load_refinement_checkpoint(
            args.resume
        )

        paths = list(
            checkpoint["input_files"]
        )

        if len(paths) != 2:
            parser.error(
                "example resume requires FOR007 "
                "and experimental input paths"
            )

        source, exp_path = map(
            Path,
            paths,
        )

        saved = checkpoint["config"].copy()
        if "refinement_bounds" in saved:
            saved["refinement_bounds"] = RefinementBounds.from_dict(saved["refinement_bounds"])
        if saved.get("bound_rules") is not None:
            saved["bound_rules"] = RefinementBoundRules.from_dict(saved["bound_rules"])

        controls = saved.pop(
            "refinement"
        ).copy()

        controls["execution"] = ForwardExecutionConfig(
            **controls["execution"]
        )

        config = FullRefinementConfig(
            LegacyRefinementConfig(
                **controls
            ),
            **saved,
        )

        state, parsed = parse_legacy_refinement(
            source,
            seed=config.refinement.seed,
            experimental_header_path=exp_path,
        )

        if config.refinement_categories is not None:
            state, _ = select_refinement_categories(
                state,
                parsed,
                config.refinement_categories,
            )

        destination = args.resume.resolve()

    # ========================================================
    # NEW RUN
    # ========================================================

    else:
        state, original = parse_legacy_refinement(
            source,
            seed=SEED,
            experimental_header_path=exp_path,
        )

        if REFINEMENT_CATEGORIES is not None:
            state, original = select_refinement_categories(
                state,
                original,
                REFINEMENT_CATEGORIES,
            )

        rules = RefinementBoundRules(SPECIFIC_BOUNDS, CATEGORY_BOUNDS,
            CATEGORY_RELATIVE_BOUNDS, CATEGORY_RELATIVE_FALLBACK)
        config = make_config(
            original,
            quick=args.quick,
            refinement_bounds=resolve_refinement_bounds(state, rules),
            bound_rules=rules,
        )

        destination = (
            ROOT
            / "tmp"
            / "full_refinement"
            / datetime.now().strftime(
                "%Y%m%d_%H%M%S_%f"
            )
        )

    # ========================================================
    # COMMAND-LINE LIMIT OVERRIDES
    # ========================================================

    if args.max_total_evaluations is not None:
        config = replace(
            config,
            max_total_evaluations=args.max_total_evaluations,
        )

    if args.max_wall_time is not None:
        config = replace(
            config,
            max_wall_time=args.max_wall_time,
        )

    # ========================================================
    # EXPERIMENTAL DATA
    # ========================================================

    records = _CaptureRecords(
        source
    )

    _parse_for007_records(
        records
    )

    row = records.values[
        "output beam count and indices"
    ]

    experimental = parse_experimental_intensities(
        exp_path,
        beam_indices_hk=tuple(
            zip(
                row[1::2],
                row[2::2],
            )
        ),
        scan=state.scan,
    )

    # ========================================================
    # PRINT CONFIGURATION
    # ========================================================

    print(
        f"Objective: {config.objective_mode}\n"
        f"Hard bounds: {config.bounds_policy}; {config.refinement_bounds}\n"
        f"Execution: {config.refinement.execution}\n"
        f"Stages: {config.refinement.stages}\n"
        f"Loops: {config.refinement.loops}\n"
        f"Seed: {config.refinement.seed}\n"
        f"Refinement categories: "
        f"{config.refinement_categories}\n"
        f"Max total evaluations: "
        f"{config.max_total_evaluations}\n"
        f"Max wall time: "
        f"{config.max_wall_time}\n"
        f"Snapshot cadence: "
        f"{config.snapshot_every}\n"
        f"Checkpoint cadence: "
        f"{config.checkpoint_every}\n"
        f"Output: {destination}",
        flush=True,
    )

    # ========================================================
    # DRY RUN
    # ========================================================

    if args.dry_run:
        layout = build_refinement_parameter_layout(
            state
        )

        print(
            f"\nActive parameters: "
            f"{len(layout.parameters)}"
        )

        for parameter in layout.parameters:
            bound = next((b for b in config.refinement_bounds.entries
                if (b.category, b.target) == (parameter.category, parameter.target)), None)
            bound_text = "unbounded" if bound is None else f"[{bound.lower}, {bound.upper}]"
            print(
                f"{parameter.vector_index + 1}: "
                f"{parameter.category} "
                f"{parameter.target}: "
                f"{parameter.initial_value} bounds={bound_text}"
            )

        print(
            "\nDry run: no objective evaluation, "
            "optimization or output directory created."
        )

        return

    # ========================================================
    # FULL REFINEMENT
    # ========================================================

    result = run_full_refinement(
        state,
        experimental,
        config,
        output_directory=destination,
        input_paths=(
            source,
            exp_path,
        ),
        resume=bool(args.resume),
    )

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    print(
        json.dumps(
            result.summary,
            indent=2,
        ),
        flush=True,
    )

    print(
        f"Output: "
        f"{result.output_directory}\n"
        f"Report: "
        f"{result.output_directory / 'report.html'}",
        flush=True,
    )


if __name__ == "__main__":
    main()

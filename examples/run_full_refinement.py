"""Edit these settings; use --quick to validate the full real-data pipeline."""
from dataclasses import replace
from datetime import datetime
from pathlib import Path
import argparse
import json

from qcbed import ForwardExecutionConfig
from qcbed.full_runner import FullRefinementConfig, run_full_refinement, load_refinement_checkpoint
from qcbed.refinement_driver import parse_legacy_refinement, LegacyRefinementConfig
from qcbed.refinement import build_refinement_parameter_layout
from qcbed.experimental import parse_experimental_intensities
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records


ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = ROOT / "reference_data" / "FOR007.DAT"
EXPERIMENTAL_FILE = ROOT / "reference_data" / "in.dat"
OBJECTIVE_MODE = "masked"
EXECUTION_BACKEND = "batched"
BATCH_SIZE = 64
OPTIMIZER_MODE = "ANSI"  # Existing ANNEAL -> AMOEBA; SPAR and other source modes also supported.
LOOPS = 1
SEED = 2026
MAXFSIMP = 20           # ANNEAL's source MAXITS, inclusive proposal limit.
MAXFPARA = 3            # PARAMIN's source cycle limit.
AMOEBA_MAX_EVALUATIONS = 20
SIMPARA_MAX_EVALUATIONS = 20
# All three controls use the SAME legacy MAXFSIMP field. For hybrid modes using
# ANNEAL and AMOEBA, their settings must agree: no per-stage override is invented.
SNAPSHOT_EVERY = 100
CHECKPOINT_EVERY = 25
MAX_TOTAL_EVALUATIONS = 1000
MAX_WALL_TIME = 3600.


def make_config(original, *, quick=False):
    cfg = replace(original, method=OPTIMIZER_MODE, loops=LOOPS, seed=SEED,
        maxfpara=MAXFPARA, execution=ForwardExecutionConfig(EXECUTION_BACKEND, BATCH_SIZE))
    limits = {"ANNEAL": MAXFSIMP, "AMOEBA": AMOEBA_MAX_EVALUATIONS, "SIMPARA": SIMPARA_MAX_EVALUATIONS}
    selected = {limits[s] for s in cfg.stages if s in limits}
    if len(selected) > 1:
        raise ValueError("selected stages share source MAXFSIMP; their configured limits must agree")
    cfg = replace(cfg, maxfsimp=next(iter(selected), MAXFSIMP))
    if quick:
        # Existing example's short ANNEAL configuration, with all physical data,
        # 18 parameters, 1521 tilts and nine beams retained. No fit retuning.
        cfg = replace(original, method="ANNE", loops=1, seed=SEED, maxfsimp=0,
            maxtemp=.001, mintemp=.0001, coolrate=.5,
            execution=ForwardExecutionConfig("batched", 64))
    return FullRefinementConfig(cfg, objective_mode="masked" if quick else OBJECTIVE_MODE,
        snapshot_every=None if quick else SNAPSHOT_EVERY,
        checkpoint_every=1 if quick else CHECKPOINT_EVERY,
        max_total_evaluations=20 if quick else MAX_TOTAL_EVALUATIONS,
        max_wall_time=600. if quick else MAX_WALL_TIME)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--quick", action="store_true")
    mode.add_argument("--resume", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-total-evaluations", type=int, help="override cumulative run budget, including on resume")
    parser.add_argument("--max-wall-time", type=float, help="override cumulative seconds, including on resume")
    args = parser.parse_args()
    source, exp_path = INPUT_FILE, EXPERIMENTAL_FILE
    if args.resume:
        checkpoint = load_refinement_checkpoint(args.resume)
        paths = list(checkpoint["input_files"])
        if len(paths) != 2:
            parser.error("example resume requires FOR007 and experimental input paths")
        source, exp_path = map(Path, paths)
        saved = checkpoint["config"].copy()
        controls = saved.pop("refinement").copy()
        controls["execution"] = ForwardExecutionConfig(**controls["execution"])
        config = FullRefinementConfig(LegacyRefinementConfig(**controls), **saved)
        state, _ = parse_legacy_refinement(source, seed=config.refinement.seed, experimental_header_path=exp_path)
        destination = args.resume.resolve()
    else:
        state, original = parse_legacy_refinement(source, seed=SEED, experimental_header_path=exp_path)
        config = make_config(original, quick=args.quick)
        destination = ROOT / "tmp" / "full_refinement" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    if args.max_total_evaluations is not None:
        config = replace(config, max_total_evaluations=args.max_total_evaluations)
    if args.max_wall_time is not None:
        config = replace(config, max_wall_time=args.max_wall_time)
    records = _CaptureRecords(source)
    _parse_for007_records(records)
    row = records.values["output beam count and indices"]
    experimental = parse_experimental_intensities(exp_path,
        beam_indices_hk=tuple(zip(row[1::2], row[2::2])), scan=state.scan)
    print(f"Objective: {config.objective_mode}\nExecution: {config.refinement.execution}\n"
        f"Stages: {config.refinement.stages}\nLoops: {config.refinement.loops}\nSeed: {config.refinement.seed}\n"
        f"Max total evaluations: {config.max_total_evaluations}\nMax wall time: {config.max_wall_time}\n"
        f"Snapshot cadence: {config.snapshot_every}\nCheckpoint cadence: {config.checkpoint_every}\n"
        f"Output: {destination}", flush=True)
    if args.dry_run:
        for p in build_refinement_parameter_layout(state).parameters:
            print(f"{p.vector_index+1}: {p.category} {p.target}: {p.initial_value}")
        print("Dry run: no objective evaluation, optimization or output directory created.")
        return
    result = run_full_refinement(state, experimental, config, output_directory=destination,
        input_paths=(source, exp_path), resume=bool(args.resume))
    print(json.dumps(result.summary, indent=2), flush=True)
    print(f"Output: {result.output_directory}\nReport: {result.output_directory / 'report.html'}", flush=True)


if __name__ == "__main__":
    main()

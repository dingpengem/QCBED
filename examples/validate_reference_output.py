"""Full-size evaluation-only reference validation; run with qcbed-python.

Writes only tmp/end_to_end_validation/. This is not historical refinement.
"""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from qcbed.experimental import parse_experimental_intensities
from qcbed.legacy_io import _parse_for007_records
from qcbed.refinement_io import _CaptureRecords
from qcbed.refinement import build_refinement_parameter_layout
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.output_io import parse_legacy_out_dat, write_legacy_out_dat
from qcbed.validation import validate_refinement
from qcbed.validation_report import generate_refinement_validation_report


def main():
    root=Path(__file__).resolve().parents[1]
    source=root/"tests/fixtures/legacy_reference"/"FOR007.DAT"
    destination=root/"tmp"/"end_to_end_validation"
    destination.mkdir(parents=True,exist_ok=True)
    state,original=parse_legacy_refinement(source,seed=2026)
    records=_CaptureRecords(source)
    _parse_for007_records(records)
    beam_record=records.values["output beam count and indices"]
    beams=tuple(zip(beam_record[1::2],beam_record[2::2]))
    experimental=parse_experimental_intensities(source.parent/state.model.experimental_intensity_filename,
                                                beam_indices_hk=beams,scan=state.scan)
    layout=build_refinement_parameter_layout(state)
    assert len(layout.parameters)==18 and experimental.intensities.shape==(1521,9)
    baseline=repr(state)
    config=replace(original,loops=0)
    started=time.perf_counter()
    print("Evaluating full 39x39 scan, nine beams, original mesh and thickness; LOOPS=0",flush=True)
    validation=validate_refinement(state,experimental,config)
    assert repr(state)==baseline
    c=validation.comparison
    write_legacy_out_dat(destination/"out.dat",c.final_simulated)
    parsed=parse_legacy_out_dat(destination/"out.dat",scan=state.scan,beam_indices_hk=beams)
    np.testing.assert_array_equal(parsed.scaled_intensity,c.final_simulated.scaled_intensity)
    np.savez_compressed(destination/"comparison.npz",experimental=experimental.intensities,
        initial_raw=c.initial_simulated.raw_intensity,final_raw=c.final_simulated.raw_intensity,
        initial_scaled=c.initial_simulated.scaled_intensity,final_scaled=c.final_simulated.scaled_intensity,
        residual=c.residual,tilts=c.final_simulated.tilts,beams=np.asarray(beams))
    summary=dict(validation_mode="evaluation-only; no optimization",active_parameters=18,
        scan_positions=1521,beams=9,
        initial_objective=float(c.initial_simulated.objective) if np.isfinite(c.initial_simulated.objective) else None,
        final_objective=float(c.final_simulated.objective) if np.isfinite(c.final_simulated.objective) else None,
        objective_status="finite" if np.isfinite(c.final_simulated.objective) else "undefined (NaN/Inf); driver cannot refine",
        zero_variance_count=int(np.count_nonzero(experimental.variances==0)),
        initial_ynorm=float(c.initial_simulated.ynorm),
        final_ynorm=float(c.final_simulated.ynorm),output_roundtrip_max_abs=0.,
        baseline_unchanged=repr(state)==baseline,
        reference_input_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    notes=["LOOPS overridden to zero: initial and final maps coincide; no refinement improvement is claimed.",
           "Diagnostic trial evaluation only. The attempted existing driver rejected the reference objective at its CHISMIN guard. "
           "Zero variances cause undefined source objective arithmetic even where weights are zero; no masking/repair was applied.",
           "tests/fixtures/legacy_reference/out.dat is empty; generated out.dat is source-derived, not a historical equivalence result."]
    # firstout.dat is a distinct pre-stage artifact, not a substitute final reference.
    historical=source.parent/"firstout.dat"
    try:
        previous=parse_legacy_out_dat(historical,scan=state.scan,beam_indices_hk=beams)
    except ValueError as exc:
        summary["firstout_comparison"]=str(exc)
    else:
        delta=c.initial_simulated.scaled_intensity-previous.scaled_intensity
        denominator=np.linalg.norm(previous.scaled_intensity)
        summary["firstout_max_abs_difference"]=float(np.max(np.abs(delta)))
        summary["firstout_relative_l2_difference"]=float(np.linalg.norm(delta)/denominator) if denominator else None
        notes.append("firstout.dat preliminary artifact comparison (not final-output equivalence): "
                     f"max absolute difference {summary['firstout_max_abs_difference']:.9g}; "
                     f"relative L2 difference {summary['firstout_relative_l2_difference']}.")
    generate_refinement_validation_report(destination/"report.html",validation,layout=layout,config=config,
                                         reference_input_name=str(source),limitations=notes)
    summary["elapsed_seconds"]=time.perf_counter()-started
    (destination/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True)
    print(destination/"report.html",flush=True)


if __name__=="__main__":
    main()

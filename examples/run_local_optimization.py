"""Edit the explicit settings below, then run with qcbed-python (no CLI)."""
from dataclasses import replace
from datetime import datetime
from pathlib import Path
import json

from qcbed.experimental import parse_experimental_intensities
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records
from qcbed.local_runner import run_local_optimization


def main():
    root=Path(__file__).resolve().parents[1]
    source=root/"tests/fixtures/legacy_reference"/"FOR007.DAT"
    state,original=parse_legacy_refinement(source,seed=2026)
    records=_CaptureRecords(source)
    _parse_for007_records(records)
    beam_record=records.values["output beam count and indices"]
    beams=tuple(zip(beam_record[1::2],beam_record[2::2]))
    experimental=parse_experimental_intensities(source.parent/state.model.experimental_intensity_filename,
                                                beam_indices_hk=beams,scan=state.scan)
    # Deliberately short validation: source MAXITS=0 still permits ONE proposal.
    # Keep all 18 parameters, 1521 tilts, nine beams and original forward physics.
    config=replace(original,method="ANNE",loops=1,maxfsimp=0,
                   maxtemp=.001,mintemp=.0001,coolrate=.5)
    destination=root/"tmp"/"local_optimization"/datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    print(f"Output: {destination}\nMasked objective; one ANNEAL proposal; seed 2026",flush=True)
    baseline=repr(state)
    run_local_optimization(state,experimental,config,objective_mode="masked",
        snapshot_every=3,output_directory=destination,reference_input_name=str(source))
    assert repr(state)==baseline
    print((destination/"summary.json").read_text(),flush=True)
    print(f"Report: {destination/'report.html'}",flush=True)


if __name__=="__main__":
    main()

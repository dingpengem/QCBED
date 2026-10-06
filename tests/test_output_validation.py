from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from qcbed.experimental import ExperimentalIntensityData, parse_experimental_intensities
from qcbed.output_io import (SimulatedCbedData, LegacySimulatedOutput, reshape_scan_map,
                             write_legacy_out_dat, parse_legacy_out_dat)
from qcbed.validation import CbedComparison, validate_refinement
from qcbed.objective import compute_legacy_objective
from qcbed.refinement import build_refinement_parameter_layout
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records
from test_refinement import state, data
from test_refinement_driver import problem, config


def comparison():
    scan=replace(state().scan,scan_columns=3,scan_rows=2,corner_b=(1.,1.),corner_c=(-1.,1.))
    x=np.arange(1,13,dtype=np.float64).reshape(6,2)
    e=ExperimentalIntensityData(3,2,0,((1,-1),(0,2)),x,np.full_like(x,2),np.full_like(x,3))
    def simulated(raw):
        objective=compute_legacy_objective(raw,e,iref=2.)
        return SimulatedCbedData(scan,e.beam_indices_hk,raw,raw,2.,0,objective.objective)
    return CbedComparison(e,simulated(x+3),simulated(x+1))


def test_explicit_row_column_reconstruction_and_roundtrip():
    c=comparison()
    for values in (c.experimental.intensities,c.initial_simulated.raw_intensity,
                   c.final_simulated.scaled_intensity,c.residual):
        for b in range(2):
            image=reshape_scan_map(values[:,b],c.final_simulated.scan)
            assert image.shape==(2,3)
            for r in range(2):
                for col in range(3):
                    assert image[r,col]==values[r*3+col,b]
            np.testing.assert_array_equal(image.ravel(),values[:,b])
            with pytest.raises(ValueError): image.setflags(write=True)
    np.testing.assert_array_equal(c.final_simulated.tilts[1],[.5,.5])
    np.testing.assert_array_equal(c.final_simulated.tilts[3],[-1,1])


def test_residual_and_beam_contributions():
    c=comparison()
    expected=c.experimental.intensities+2
    np.testing.assert_array_equal(c.residual,expected)
    for b in range(2):
        assert c.numerator_by_beam[b]==sum(3*x*x/2 for x in expected[:,b])
        assert c.denominator_by_beam[b]==sum(3*x*x/2 for x in c.experimental.intensities[:,b])
    total=np.sqrt(sum(c.numerator_by_beam)/sum(c.denominator_by_beam))
    assert total==pytest.approx(c.final_simulated.objective,rel=1e-15)
    for a in (c.residual,c.final_simulated.raw_intensity,c.final_simulated.objective_intensity,
              c.final_simulated.scaled_intensity,c.numerator_by_beam):
        with pytest.raises(ValueError): a.setflags(write=True)


def test_output_roundtrip_order_and_wrapping(tmp_path):
    c=comparison()
    path=tmp_path/"out.dat"
    write_legacy_out_dat(path,c.final_simulated)
    lines=path.read_text().splitlines()
    assert len(lines)==8 and lines[1]==""
    assert lines[0]=="4.00000000000000000E+00 8.00000000000000000E+00 1.20000000000000000E+01"
    # Runtime-dependent wrapping may split a logical row over physical records.
    path.write_text("\n".join(path.read_text().replace("E","D").split()))
    parsed=parse_legacy_out_dat(path,scan=c.final_simulated.scan,beam_indices_hk=c.experimental.beam_indices_hk)
    assert parsed.beam_indices_hk==((1,-1),(0,2))
    np.testing.assert_array_equal(parsed.scaled_intensity,c.final_simulated.scaled_intensity)


@pytest.mark.parametrize("text",["", "1", " ".join(["1"]*13), " ".join(["NaN"]*12)])
def test_invalid_output_rejected(tmp_path,text):
    path=tmp_path/"invalid.dat"
    path.write_text(text)
    c=comparison()
    with pytest.raises(ValueError,match="invalid.dat: section 'simulated output'"):
        parse_legacy_out_dat(path,scan=c.final_simulated.scan,beam_indices_hk=c.experimental.beam_indices_hk)


def test_gradient_output_uses_objective_values_not_raw():
    c=comparison()
    raw=np.full((6,2),10.)
    processed=np.full((6,2),-3.)
    s=SimulatedCbedData(c.final_simulated.scan,c.experimental.beam_indices_hk,raw,processed,2.,1,.5)
    np.testing.assert_array_equal(s.scaled_intensity,np.full((6,2),-6.))
    assert s.raw_intensity[0,0]==10


def test_reference_dimensions_and_maps():
    source=Path(__file__).parents[1]/"tests/fixtures/legacy_reference"/"FOR007.DAT"
    s,config= parse_legacy_refinement(source,seed=0)
    records=_CaptureRecords(source)
    _parse_for007_records(records)
    row=records.values["output beam count and indices"]
    beams=tuple(zip(row[1::2],row[2::2]))
    e=parse_experimental_intensities(source.parent/"in.dat",beam_indices_hk=beams,scan=s.scan)
    assert e.intensities.shape==(1521,9) and len(build_refinement_parameter_layout(s).parameters)==18
    for b in range(9):
        image=reshape_scan_map(e.intensities[:,b],s.scan)
        assert image.shape==(39,39) and image[7,19]==e.intensities[7*39+19,b]
    assert (source.parent/"out.dat").stat().st_size==0
    assert np.count_nonzero(e.variances==0)==117
    # Source 0/0 is not repaired even if the corresponding weight is zero.
    assert np.isnan(compute_legacy_objective(np.ones((1521,9)),e,iref=s.iref).objective)


def test_real_short_refinement_and_evaluation_only():
    s=problem()
    original=repr(s)
    for loops in (0,1):
        validation=validate_refinement(s,data(),config("ANSI",loops=loops,maxfsimp=0))
        r,c=validation.refinement,validation.comparison
        assert c.initial_simulated.objective==r.initial_objective
        assert c.final_simulated.objective==pytest.approx(r.final_objective,rel=1e-12)
        assert repr(s)==original
        if loops==0:
            np.testing.assert_array_equal(c.initial_simulated.raw_intensity,c.final_simulated.raw_intensity)


def test_metadata_mismatch_rejected():
    c=comparison()
    wrong=replace(c.final_simulated,beam_indices_hk=((0,2),(1,-1)))
    with pytest.raises(ValueError,match="beam order"):
        CbedComparison(c.experimental,c.initial_simulated,wrong)


def test_evaluation_only_preserves_undefined_objective():
    s=problem()
    e=data()
    e=replace(e,variances=np.zeros_like(e.variances),weights=np.zeros_like(e.weights))
    # Fixed finite YNORM leaves maps valid while the objective is undefined.
    s=replace(s,iref=1.)
    v=validate_refinement(s,e,config("ANSI",loops=0))
    assert np.isnan(v.comparison.final_simulated.objective)
    assert np.all(np.isfinite(v.comparison.final_simulated.scaled_intensity))
    assert np.all(np.isnan(v.comparison.numerator_by_beam))
    assert v.refinement.evaluations==0 and v.refinement.stage_history==()
    with pytest.raises(ValueError,match="CHISMIN"):
        validate_refinement(s,e,config("ANSI",loops=1))

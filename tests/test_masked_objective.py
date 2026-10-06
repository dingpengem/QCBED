from pathlib import Path

import numpy as np
import pytest

from qcbed.masked_objective import compute_masked_objective, compute_selected_objective
from qcbed.objective import compute_legacy_objective
from test_objective import data


@pytest.mark.parametrize("iref,scale,error", [(1.,1.,1.5),(0.,-6/13,np.sqrt(56/65)),
                                            (-1.,10/13,np.sqrt(24/13))])
def test_independent_weighted_hand_case(iref,scale,error):
    c=np.array([[1.,2.],[3.,4.]])
    # Common fractional weight factor preserves hand-calculated scale/error,
    # halves the normalizer. The masked variance is deliberately zero.
    e=data([[2,-4],[6,0]],[[1,2],[0,4]],[[.5,1],[0,1]])
    before=[a.copy() for a in (c,e.intensities,e.variances,e.weights)]
    with np.errstate(all="raise"):
        r=compute_masked_objective(c,e,iref=iref)
    assert r.scaling==pytest.approx(scale,rel=1e-15)
    assert r.objective==pytest.approx(error,rel=1e-15)
    assert r.weighted_experimental_norm==10
    if iref==1: assert r.weighted_squared_error==22.5
    assert np.isnan(compute_selected_objective(c,e,iref=iref).objective)
    assert compute_selected_objective(c,e,iref=iref,objective_mode="masked")==r
    for a,b in zip((c,e.intensities,e.variances,e.weights),before):
        np.testing.assert_array_equal(a,b)


@pytest.mark.parametrize("variance", [0.,-1.])
@pytest.mark.parametrize("iref", [1.,0.,-1.])
def test_visible_invalid_variance_always_rejected(variance,iref):
    e=data([[2]],[[variance]],[[.25]])
    with pytest.raises(ValueError,match="nonzero-weight.*finite positive variance.*tilt 0, beam 0"):
        compute_masked_objective(np.ones((1,1)),e,iref=iref)


@pytest.mark.parametrize("iref", [1.,0.,-1.])
def test_matches_legacy_when_all_variances_positive(iref):
    c=np.array([[1.,2.],[3.,4.]])
    e=data([[2,-4],[6,0]],[[1,2],[3,4]],[[.5,1],[0,1]])
    a=compute_masked_objective(c,e,iref=iref)
    b=compute_legacy_objective(c,e,iref=iref)
    for field in a.__dataclass_fields__: assert getattr(a,field)==getattr(b,field)


def test_empty_mask_and_singular_scaling_and_mode():
    c=np.ones((1,1))
    with pytest.raises(ValueError,match="at least one"):
        compute_masked_objective(c,data([[1]],[[0]],[[0]]),iref=1)
    with pytest.raises(ValueError,match="YNORM denominator"):
        compute_masked_objective(c*0,data([[1]]),iref=0)
    with pytest.raises(ValueError,match="normalizer"):
        compute_masked_objective(c,data([[0]]),iref=1)
    with pytest.raises(ValueError,match="objective_mode"):
        compute_selected_objective(c,data([[1]]),iref=1,objective_mode="auto")


def test_reference_mask_gives_finite_score_without_changing_legacy():
    from qcbed.experimental import parse_experimental_intensities
    # The objective does not depend on H/K identities, only preserved beam order.
    p=Path(__file__).parents[1]/"tests/fixtures/legacy_reference"/"in.dat"
    beams=((0,0),(1,1),(1,-1),(-1,-1),(-1,1),(2,0),(0,-2),(-2,0),(0,2))
    e=parse_experimental_intensities(p,beam_indices_hk=beams)
    assert e.intensities.shape==(1521,9)
    c=np.ones_like(e.intensities)
    for iref in (1.,0.,-1.):
        assert np.isfinite(compute_masked_objective(c,e,iref=iref).objective)
        assert np.isnan(compute_legacy_objective(c,e,iref=iref).objective)

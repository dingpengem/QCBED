import numpy as np
import pytest
from qcbed.depth_geometry import requested_depth_to_sample_geometry as geometry


def test_source_loop_and_boundary_plateaus():
    baseline=geometry(840.,(2.0245,2.0245))
    assert baseline.occurrence_count==415
    assert baseline.effective_depth==840.1674999999964
    assert baseline.stored_slice_count==2
    lo,hi=baseline.plateau_lower_exclusive,baseline.effective_depth
    assert geometry(lo,(2.0245,2.0245)).occurrence_count==414
    assert geometry(np.nextafter(lo,np.inf),(2.0245,2.0245)).sequence==baseline.sequence
    assert geometry(hi,(2.0245,2.0245)).sequence==baseline.sequence
    assert geometry(np.nextafter(hi,np.inf),(2.0245,2.0245)).occurrence_count==416
    assert geometry(839.,(2.0245,2.0245)).sequence==baseline.sequence


@pytest.mark.parametrize('thicknesses',[(1.,),(1.,2.,.5),(.1,.2)])
def test_mixed_patterns_preserve_literal_sum(thicknesses):
    for request in (.01,.3,1.,3.,7.77):
        sequence=[];total=0.
        while total<request:
            i=len(sequence)%len(thicknesses);sequence.append(i);total+=thicknesses[i]
        observed=geometry(request,thicknesses)
        assert observed.sequence==tuple(sequence) and observed.effective_depth==total
        assert observed.plateau_lower_exclusive<request<=observed.effective_depth


@pytest.mark.parametrize('depth,thicknesses',[(0,(1,)),(-1,(1,)),(np.nan,(1,)),(np.inf,(1,)),
    (1,()),(1,(0,)),(1,(-1,)),(1,(np.inf,)),(True,(1,))])
def test_invalid_depth(depth,thicknesses):
    with pytest.raises(ValueError):geometry(depth,thicknesses)

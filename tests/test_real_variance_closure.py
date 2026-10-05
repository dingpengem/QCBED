"""Real sigma replay with archived A7/A8/A10 geometry; no physical HK inference."""
import numpy as np
from examples.a9_real_variance_closure import prepare_reference_variance


def test_real_nine_beam_variance_closure():
    datasets,summary=prepare_reference_variance()
    expected=((0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))
    assert summary['sigma_shape']==[256,256] and summary['sigma_values']==65536
    for name,count in [('conventional',886),('gradient',809)]:
        data=datasets[name];stats=summary['modes'][name]['statistics']
        assert data.variances.shape==data.intensities.shape==data.weights.shape==(1444,9)
        assert data.lattice_coordinates==expected
        assert stats['support_per_beam']==[count]*9
        assert stats['active_count']==count*9
        assert stats['active_min']>0
        assert stats['active_nonpositive_count']==stats['nan_count']==stats['inf_count']==0
        assert stats['finite_count']==12996
        assert data.scan_geometry==datasets['conventional'].scan_geometry
        assert np.all(np.isfinite(data.variances))
        assert all(b.indices is None for b in data.beams) # physical HK remains unresolved
        assert data.variances.flags.writeable is False
    assert summary['a4_refit'] is False and summary['objective']=='not evaluated'

import numpy as np
import pytest
from qcbed.masked_pattern_diagnostics import masked_comparison,save_masked_pattern_match


def test_exact_values_mask_orientation_and_plot(tmp_path):
    x=np.arange(24,dtype=float).reshape(6,4)
    w=np.ones_like(x);w[1,2]=0
    arrays=dict(experimental=x,scaled_calculated=x*2,residual=x,weights=w)
    before={k:v.copy() for k,v in arrays.items()}
    data=save_masked_pattern_match(arrays,2,3,['a','b','a','d'],tmp_path/'plot.png',title='test',ynorm=2.)
    assert data['experimental'].shape==(4,2,3)
    assert np.isnan(data['experimental'][2,0,1])
    mask=data['mask'];expected=(x*2).T.reshape(4,2,3)
    np.testing.assert_array_equal(data['scaled_calculated'][mask],expected[mask])
    for k,v in before.items():np.testing.assert_array_equal(arrays[k],v)
    saved=np.load(tmp_path/'plot.npz')
    np.testing.assert_array_equal(saved['mask'],mask)
    np.testing.assert_array_equal(saved['residual'],data['residual'])
    assert saved['ynorm']==2 and (tmp_path/'plot.png').stat().st_size>0
    with pytest.raises(ValueError):masked_comparison(arrays,3,3)

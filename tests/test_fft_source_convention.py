"""Independent source-sign oracles; no historical executable claim."""
from pathlib import Path

import numpy as np
import pytest

from qcbed.wave_propagation import real_to_reciprocal_wave, reciprocal_to_real_wave
from qcbed.fft_backend import batched_wave_real_to_reciprocal, batched_wave_reciprocal_to_real
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement_trial import _prepare_trial
from qcbed.simulation import SamplePreparationParameters, prepare_multislice_sample, build_tilt_propagators, simulate_one_tilt
from test_multislice import nested_transform, nested_sequence


def test_four_point_literal_butterfly_and_orientation():
    # fftsg.f:2858-2876 CFTB040: [0,1,0,0] -> [1,-i,-1,+i].
    wave = np.array([[0,1,0,0]],dtype=np.complex128)
    expected = np.array([[.25,-.25j,-.25,.25j]])
    np.testing.assert_array_equal(real_to_reciprocal_wave(wave,ordering='unshifted'), expected)
    np.testing.assert_array_equal(reciprocal_to_real_wave(expected,ordering='unshifted'),wave)
    # A complex asymmetric input distinguishes labeled +H/-H intensities.
    wave = np.array([[1,1j,0,0]],dtype=np.complex128)
    actual = real_to_reciprocal_wave(wave,ordering='unshifted')
    np.testing.assert_array_equal(actual, [[.25+.25j,.5,.25-.25j,0]])
    np.testing.assert_array_equal(actual.real**2+actual.imag**2,[[.125,.25,.125,0]])


@pytest.mark.parametrize('shape',[(3,5),(4,8),(1,4)])
def test_each_direction_and_batch_against_scalar_sums(shape):
    n=np.prod(shape)
    a=(np.arange(n).reshape(shape)/13+1j*np.arange(n)[::-1].reshape(shape)/17).astype(np.complex128)
    for to_reciprocal,serial,batch in (
        (True,real_to_reciprocal_wave,batched_wave_real_to_reciprocal),
        (False,reciprocal_to_real_wave,batched_wave_reciprocal_to_real)):
        expected=nested_transform(a,to_reciprocal=to_reciprocal)
        np.testing.assert_allclose(serial(a,ordering='unshifted'),expected,atol=7e-14,rtol=0)
        np.testing.assert_allclose(batch(np.stack([a,a*1j,a]))[0],expected,atol=7e-14,rtol=0)
    back=reciprocal_to_real_wave(real_to_reciprocal_wave(a,ordering='unshifted'),ordering='unshifted')
    np.testing.assert_allclose(back,a,atol=2e-15,rtol=0)


@pytest.mark.parametrize('h,k',[(0.,0.),(.2,-.3),(-.2,.3)])
def test_real_al_input_three_slices_against_nested_dft(h,k):
    # Real geometry/species/atoms/absorption and 16x16 mesh; short sequence
    # keeps the independent O(N^2) scalar DFT practical in the unit suite.
    s,_=parse_legacy_refinement(Path(__file__).parents[1]/'reference_data'/'FOR007.DAT',seed=2026)
    applied,_,_,coefficients=_prepare_trial(s,original_voltage=s.model.energy_kv)
    beams=((0,0),(1,1),(1,-1),(-1,-1),(-1,1),(2,0),(0,-2),(-2,0),(0,2))
    p=SamplePreparationParameters(s.model.energy_kv,applied,s.dsap,coefficients,(0,)*3,
        beams,'unshifted','fft','legacy','unity')
    prepared=prepare_multislice_sample(s.model,p)
    kernels=build_tilt_propagators(prepared,tilt_h_undilated_fractional=h,tilt_k_undilated_fractional=k)
    expected=nested_sequence(np.ones((16,16),dtype=np.complex128),
        [a.phase_grating for a in prepared.slice_types],kernels,(0,)*3)
    result=simulate_one_tilt(prepared,tilt_h_undilated_fractional=h,tilt_k_undilated_fractional=k)
    amplitudes=np.array([expected[bk%16,bh%16] for bh,bk in beams])
    np.testing.assert_allclose(result.final_reciprocal_wave,expected,atol=8e-15,rtol=0)
    np.testing.assert_allclose(result.amplitudes,amplitudes,atol=8e-15,rtol=0)
    np.testing.assert_allclose(result.intensities,abs(amplitudes)**2,atol=1e-14,rtol=0)

import cmath
from dataclasses import FrozenInstanceError, replace
import math
from pathlib import Path
import weakref

import numpy as np
import pytest

from qcbed import (
    AtomSite, Species, MeshConfig, SliceGeometry, AbsorptionCoefficients,
    SliceAppliedParameters, SamplePreparationParameters, prepare_multislice_sample,
    simulate_one_tilt, TiltSeriesResult, simulate_tilt_series_serial,
    simulate_tilt_series_batched, parse_for007,
)
from qcbed.legacy_io import LegacyInput
from qcbed.fft_backend import (
    batched_real_to_reciprocal_coefficients, batched_reciprocal_to_real_sum,
)


TILTS = np.array([[0., 0.], [.2, 0.], [-.3, 0.], [0., .4], [0., -.5], [.2, -.1], [.2, 0.]])


def sample(*, mesh=(6, 4), ordering="unshifted", method="fft", absorption=False,
           dilation=True, gamma=60, dc=False):
    geometries = (SliceGeometry(2, .125, -.25, 20, 10, gamma),
                  SliceGeometry(3, -.3, .2, 20, 10, gamma))
    model = LegacyInput("tilts", 200, False, (0, 0, 1), MeshConfig(*mesh), geometries,
        (Species("Al", 13),), ((AtomSite(1, .25, .125, .5, "S"),),
                             (AtomSite(1, .375, .25, .75, "S"),)),
        "out.dat", "in.dat", title_raw="tilts", absorption_raw="YES" if absorption else "NO")
    parameters = SamplePreparationParameters(200,
        (SliceAppliedParameters((.86,), .1, 1.1 if dilation else 1),
         SliceAppliedParameters((1.2,), .2, .95 if dilation else 1)),
        0 if dc else .5, (AbsorptionCoefficients(0, 1, 1),) if absorption else None,
        (0, 1, 0, 1), ((1, 0), (0, 0), (-1, 0), (1, 0), (0, -1)),
        ordering, method, "legacy", "unity")
    return prepare_multislice_sample(model, parameters)


@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
@pytest.mark.parametrize("mesh", [(4, 4), (6, 4)])
@pytest.mark.parametrize("method", ["direct", "fft"])
@pytest.mark.parametrize("absorption", [False, True])
@pytest.mark.parametrize("batch_size", [1, 2, 20])
def test_serial_batched_matrix(ordering, mesh, method, absorption, batch_size):
    prepared = sample(mesh=mesh, ordering=ordering, method=method, absorption=absorption)
    before = tuple(s.phase_grating.copy() for s in prepared.slice_types)
    serial = simulate_tilt_series_serial(prepared, TILTS)
    batched = simulate_tilt_series_batched(prepared, TILTS, batch_size=batch_size)
    repeated = np.stack([simulate_one_tilt(prepared, tilt_h_undilated_fractional=h,
                        tilt_k_undilated_fractional=k).amplitudes for h, k in TILTS])
    np.testing.assert_array_equal(serial.amplitudes, repeated)
    np.testing.assert_allclose(batched.amplitudes, serial.amplitudes, atol=3e-15, rtol=0)
    np.testing.assert_allclose(batched.intensities, serial.intensities, atol=4e-15, rtol=0)
    for result in (serial, batched):
        np.testing.assert_array_equal(result.tilts, TILTS)
        np.testing.assert_array_equal(result.amplitudes[1], result.amplitudes[-1])
        np.testing.assert_array_equal(result.amplitudes[:, 0], result.amplitudes[:, 3])
        assert result.beam_indices_hk == prepared.parameters.beam_indices_hk
        assert result.amplitudes.dtype == np.complex128
        assert result.intensities.dtype == result.tilts.dtype == np.float64
        assert result.amplitudes.shape == (7, 5)
    for stored, old in zip(prepared.slice_types, before):
        np.testing.assert_array_equal(stored.phase_grating, old)
    assert serial.execution_mode == "serial" and serial.batch_size is None
    assert batched.execution_mode == "batched" and batched.batch_size == batch_size


@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_no_cross_talk(ordering):
    prepared = sample(mesh=(6, 4), ordering=ordering)
    tilts = np.array([[.1, -.2], [.3, .4], [-.5, .6]])
    before = simulate_tilt_series_batched(prepared, tilts, batch_size=3)
    tilts[1] = [7., -9.]
    after = simulate_tilt_series_batched(prepared, tilts, batch_size=3)
    np.testing.assert_array_equal(before.amplitudes[[0, 2]], after.amplitudes[[0, 2]])
    np.testing.assert_array_equal(before.intensities[[0, 2]], after.intensities[[0, 2]])
    assert np.max(abs(before.amplitudes[1]-after.amplitudes[1])) > 1e-8
    np.testing.assert_array_equal(before.tilts[1], [.3, .4])


def dft(values, to_reciprocal):
    ny, nx = values.shape
    return np.array([[sum(values[j, i]*cmath.exp((-1 if to_reciprocal else 1)*2j*math.pi*
        (i*c/nx+j*r/ny)) for j in range(ny) for i in range(nx))/(nx*ny if to_reciprocal else 1)
        for c in range(nx)] for r in range(ny)], dtype=np.complex128)


def independent_tiny(h0, k0, *, return_states=False):
    """Printed Al coefficients, Cartesian dual vectors, scalar nested DFTs.

    2x4 mask retains K=0,H=-1,0,1 only: STL=0,.025. No production
    grating, propagator, scattering, grid or FFT is used for expectations.
    """
    gratings, kernels = [], []
    for x, occupancy, vacancy, b, z, shift in ((.25, .5, .1, .86, 2, .125),
                                             (.375, .75, .2, 1.2, 3, -.3)):
        coefficients = {h: (281.953595001 if h == 0 else 275.340445007875)/200
            *occupancy*(1-vacancy)*math.exp(-b*(h*.025)**2)*cmath.exp(1j*6.283185307*h*x)
            for h in (-1, 0, 1)}
        grating = np.array([[cmath.exp(-1j*.0007288360746343612*sum(
            f*cmath.exp(-2j*math.pi*h*c/4) for h, f in coefficients.items()))
            for c in range(4)] for _ in range(2)])
        kernel = np.zeros((2, 4), dtype=np.complex128)
        for h in (-1, 0, 1):
            gx, gy = h/20, -h*math.cos(3.14159265358/2)/20
            tx, ty = h0/20, -h0*math.cos(3.14159265358/2)/20+k0/10
            phase = 6.2831853*(.0250794525129636/2*(gx*gx+gy*gy+2*(tx*gx+ty*gy))*z-h*shift)
            kernel[0, h % 4] = cmath.exp(1j*phase)
        gratings.append(grating)
        kernels.append(kernel)
    wave = np.ones((2, 4), dtype=np.complex128)
    transmitted_states, reciprocal_states = [], []
    for occurrence, index in enumerate((0, 1, 0, 1)):
        transmitted = wave*gratings[index]
        reciprocal = dft(transmitted, True)*kernels[index]
        transmitted_states.append(transmitted)
        reciprocal_states.append(reciprocal)
        if occurrence < 3:
            wave = dft(reciprocal, False)
    if return_states:
        return transmitted_states, reciprocal_states
    return np.array([reciprocal[0, 1], reciprocal[0, 0], reciprocal[0, 3],
                     reciprocal[0, 1], reciprocal[1, 0]])


@pytest.mark.parametrize("method", ["direct", "fft"])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_independent_nested_dft(method, ordering):
    prepared = sample(mesh=(4, 2), gamma=90, dilation=False, method=method, ordering=ordering)
    expected = np.stack([independent_tiny(h, k) for h, k in TILTS])
    for result in (simulate_tilt_series_serial(prepared, TILTS),
                   simulate_tilt_series_batched(prepared, TILTS, batch_size=3)):
        np.testing.assert_allclose(result.amplitudes, expected, atol=3e-15, rtol=0)
        np.testing.assert_allclose(result.intensities, abs(expected)**2, atol=4e-15, rtol=0)


def test_intermediate_complex_waves_against_nested_dft(monkeypatch):
    import qcbed.tilt_series as module
    prepared = sample(mesh=(4, 2), gamma=90, dilation=False, ordering="centered")
    tilts = TILTS[[0, 3, 5]]
    expected = [independent_tiny(h, k, return_states=True) for h, k in tilts]
    calls = [0, 0]
    for name, which in (("batched_wave_real_to_reciprocal", 0),
                         ("batched_wave_reciprocal_to_real", 1)):
        original = getattr(module, name)
        def transform(values, _original=original, _which=which):
            reference = np.stack([states[_which][calls[_which]] for states in expected])
            np.testing.assert_allclose(values, reference, atol=3e-15, rtol=0)
            calls[_which] += 1
            return _original(values)
        monkeypatch.setattr(module, name, transform)
    simulate_tilt_series_batched(prepared, tilts, batch_size=3)
    assert calls == [4, 3]


def test_one_zero_tilt_no_extra_rmt_and_immutable_result():
    prepared = sample(dc=True, absorption=False)
    tilts = np.zeros((1, 2), dtype=np.float64)
    result = simulate_tilt_series_batched(prepared, tilts, batch_size=np.int64(8))
    # A spatially uniform, elastic-only grating changes phase, not intensity.
    np.testing.assert_allclose(result.intensities, [[0, 1, 0, 0, 0]], atol=2e-15, rtol=0)
    tilts[:] = 9
    np.testing.assert_array_equal(result.tilts, [[0, 0]])
    for value in (result.tilts, result.amplitudes, result.intensities):
        with pytest.raises(ValueError):
            value.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        result.batch_size = 1


@pytest.mark.parametrize("bad", [[], np.zeros((0, 2)), np.zeros(2), np.zeros((2, 3)),
    np.zeros((2, 1, 2)), np.array([[np.nan, 0]]), np.array([[0, np.inf]]),
    np.array([[-np.inf, 0]]), np.zeros((2, 2), dtype=int), np.zeros((2, 2), dtype=bool),
    np.zeros((2, 2), dtype=np.float32), np.zeros((2, 2), dtype=np.complex128)])
def test_invalid_tilts(bad):
    prepared = sample()
    with pytest.raises((TypeError, ValueError), match="tilts"):
        simulate_tilt_series_serial(prepared, bad)
    with pytest.raises((TypeError, ValueError), match="tilts"):
        simulate_tilt_series_batched(prepared, bad, batch_size=2)


@pytest.mark.parametrize("bad", [True, np.bool_(False), 0, -1, 1., .5, np.nan, np.inf, "2", None])
def test_invalid_batch_size(bad):
    with pytest.raises((TypeError, ValueError), match="batch_size"):
        simulate_tilt_series_batched(sample(), TILTS, batch_size=bad)


def test_bad_prepared_and_result_validation():
    with pytest.raises(TypeError, match="prepared_sample"):
        simulate_tilt_series_serial(None, TILTS)
    with pytest.raises(TypeError, match="prepared_sample"):
        simulate_tilt_series_batched(None, TILTS, batch_size=2)
    result = simulate_tilt_series_serial(sample(), TILTS)
    for changes, match in ((dict(amplitudes=np.ones((7, 5))), "complex128"),
        (dict(amplitudes=np.ones((1, 5), dtype=np.complex128)), "shape"),
        (dict(amplitudes=np.full((7, 5), np.inf, dtype=np.complex128)), "finite"),
        (dict(execution_mode="auto"), "execution_mode"), (dict(batch_size=1), "serial")):
        with pytest.raises((TypeError, ValueError), match=match):
            replace(result, **changes)


def test_empty_beams_and_noncontiguous_input():
    prepared = sample()
    prepared = replace(prepared, parameters=replace(prepared.parameters, beam_indices_hk=()))
    tilts = TILTS[::-2]
    for result in (simulate_tilt_series_serial(prepared, tilts),
                   simulate_tilt_series_batched(prepared, tilts, batch_size=2)):
        assert result.amplitudes.shape == result.intensities.shape == (4, 0)
        np.testing.assert_array_equal(result.tilts, tilts)


def test_chunk_lifetime_caching_and_no_final_inverse(monkeypatch):
    import qcbed.tilt_series as module
    prepared = sample(ordering="centered")
    counts = dict(mask=0, kernel=0, forward=0, inverse=0)
    references = []
    for name, key in (("propagation_boundary_mask", "mask"), ("propagation_kernel", "kernel"),
        ("batched_wave_real_to_reciprocal", "forward"),
        ("batched_wave_reciprocal_to_real", "inverse")):
        original = getattr(module, name)
        def wrapper(*args, _original=original, _key=key, **kwargs):
            counts[_key] += 1
            if _key in ("forward", "inverse"):
                assert args[0].shape[0] <= 2
                references.append(weakref.ref(args[0]))
            if _key == "kernel":
                assert kwargs["boundary"] == "none"
            return _original(*args, **kwargs)
        monkeypatch.setattr(module, name, wrapper)
    original_chunk = module._chunk
    def chunk(*args):
        assert all(ref() is None for ref in references)
        return original_chunk(*args)
    monkeypatch.setattr(module, "_chunk", chunk)
    simulate_tilt_series_batched(prepared, TILTS, batch_size=2)
    assert counts == dict(mask=1, kernel=14, forward=16, inverse=12)
    assert all(ref() is None for ref in references)


def test_batched_backend_against_independent_dft():
    values = (np.arange(36).reshape(3, 3, 4)/13 + 1j*np.arange(36)[::-1].reshape(3, 3, 4)/17).astype(np.complex128)
    for function, positive in ((batched_real_to_reciprocal_coefficients, True),
                                (batched_reciprocal_to_real_sum, False)):
        expected = np.stack([dft(v.conj(), positive).conj() for v in values])
        np.testing.assert_allclose(function(values), expected, atol=4e-14, rtol=0)
    np.testing.assert_allclose(batched_reciprocal_to_real_sum(
        batched_real_to_reciprocal_coefficients(values)), values, atol=1e-15, rtol=0)


@pytest.mark.parametrize("bad", [np.zeros((2, 3)), np.zeros((2, 3, 4)),
    np.zeros((0, 3, 4), dtype=np.complex128), np.full((2, 3, 4), np.nan, dtype=np.complex128)])
def test_batch_backend_validation(bad):
    for function in (batched_real_to_reciprocal_coefficients, batched_reciprocal_to_real_sum):
        with pytest.raises((TypeError, ValueError)):
            function(bad)


def reference_sample():
    path = Path(__file__).parents[1]/"tests/fixtures/legacy_reference"/"FOR007.DAT"
    model = parse_for007(path)
    lines = path.read_text().splitlines()
    assert list(map(float, lines[10].split())) == [840., 1.]
    assert float(lines[17]) == .86 and float(lines[7]) == 2.5
    assert list(map(float, lines[13].split())) == [0., 0.]
    assert list(map(float, lines[14].split())) == [0., 0.]
    assert list(map(float, lines[16].split())) == [0., 0., 1., 1., 1., 1.]
    assert all(list(map(float, line.split()[2:])) == [1., 1.] for line in lines[46:54])
    assert 414*model.slices[0].thickness_angstrom < 840 <= 415*model.slices[0].thickness_angstrom
    start = next(i for i, line in enumerate(lines) if line.strip() == "out.dat")+2
    beams = tuple(tuple(map(int, line.split())) for line in lines[start:start+9])
    parameters = SamplePreparationParameters(200, (SliceAppliedParameters((.86,), 0, 1),),
        2.5, (AbsorptionCoefficients(0, 1, 1),), (0,)*415, beams, "unshifted", "fft", "legacy", "unity")
    return prepare_multislice_sample(model, parameters)


def test_real_al_subset():
    prepared = reference_sample()
    tilts = np.array([[0., 0.], [.1, -.2], [-.2, .15], [.1, -.2], [.3, .2]])
    serial = simulate_tilt_series_serial(prepared, tilts)
    batched = simulate_tilt_series_batched(prepared, tilts, batch_size=2)
    np.testing.assert_allclose(batched.amplitudes, serial.amplitudes, atol=3e-14, rtol=0)
    np.testing.assert_allclose(batched.intensities, serial.intensities, atol=3e-14, rtol=0)
    assert batched.amplitudes.shape == (5, 9)
    assert np.all(np.isfinite(batched.amplitudes)) and np.all(np.isfinite(batched.intensities))
    assert np.all(batched.intensities >= 0)
    np.testing.assert_array_equal(batched.tilts, tilts)
    assert batched.beam_indices_hk == prepared.parameters.beam_indices_hk

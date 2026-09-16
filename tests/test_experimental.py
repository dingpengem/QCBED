from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import numpy as np
import pytest

from qcbed import (ExperimentalIntensityData, parse_experimental_intensities,
                   LegacyScanDefinition, parse_for007_scan, generate_scan_tilts)
from qcbed.legacy_io import LegacyInputParseError


BEAMS = ((2, -1), (-3, 0), (2, -1))
REFERENCE_BEAMS = ((0, 0), (1, 1), (1, -1), (-1, -1), (-1, 1),
                   (2, 0), (0, -2), (-2, 0), (0, 2))
REFERENCE = Path(__file__).parents[1]/"reference_data"


def write_data(tmp_path):
    path = tmp_path/"experimental.dat"
    path.write_text("2 2 1\n-1 0\n3 4\n10 20\n30 40\n5 6\n7 8\n"
                    "0 2\n3 4\n5 6\n7 8\n9 10\n11 12\n"
                    "0 1\n.5 1\n1 0\n1 .25\n1 1\n0 0\n")
    return path


def test_three_blocks_orientation_values_and_metadata(tmp_path):
    data = parse_experimental_intensities(write_data(tmp_path), beam_indices_hk=BEAMS)
    np.testing.assert_array_equal(data.intensities, [[-1, 10, 5], [0, 20, 6], [3, 30, 7], [4, 40, 8]])
    np.testing.assert_array_equal(data.variances, [[0, 5, 9], [2, 6, 10], [3, 7, 11], [4, 8, 12]])
    np.testing.assert_array_equal(data.weights, [[0, 1, 1], [1, 0, 1], [.5, 1, 0], [1, .25, 0]])
    assert data.beam_indices_hk == BEAMS and data.gradient_option == 1
    assert data.scan_columns == data.scan_rows == 2
    for array in (data.intensities, data.variances, data.weights):
        assert array.dtype == np.float64 and array.shape == (4, 3)
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        data.scan_columns = 3


def test_fortran_syntax_continuation_and_advancing_read(tmp_path):
    path = tmp_path/"syntax.dat"
    path.write_text("2,\n1,0,999\n+1D1,\n-2d-1,777\n2*0\n-1,2.5\n\n")
    data = parse_experimental_intensities(path, beam_indices_hk=((0, 0),))
    np.testing.assert_array_equal(data.intensities[:, 0], [10, -.2])
    np.testing.assert_array_equal(data.variances[:, 0], [0, 0])
    np.testing.assert_array_equal(data.weights[:, 0], [-1, 2.5])
    # Signed weights and zero variances are retained, not validated for fitting.


@pytest.mark.parametrize("token", ["NaN", "Inf", "-Inf", "1+", "word", "1/", "'2'"])
@pytest.mark.parametrize("line", [1, 7, 13])
def test_invalid_numeric_values_with_context(tmp_path, token, line):
    path = write_data(tmp_path)
    records = path.read_text().splitlines()
    records[line] = f"{token} 0"
    path.write_text("\n".join(records))
    with pytest.raises(LegacyInputParseError) as error:
        parse_experimental_intensities(path, beam_indices_hk=BEAMS)
    message = str(error.value)
    assert str(path) in message and "section 'experimental" in message
    assert f"record {line+1}" in message and "expected" in message


@pytest.mark.parametrize("header", ["0 2 0", "2 -1 0", "2.0 2 0", "2 2 NaN", "2 2 .5"])
def test_invalid_header(tmp_path, header):
    path = write_data(tmp_path)
    path.write_text(header+"\n"+"\n".join(path.read_text().splitlines()[1:]))
    with pytest.raises(LegacyInputParseError, match="experimental header.*record 1"):
        parse_experimental_intensities(path, beam_indices_hk=BEAMS)


@pytest.mark.parametrize("change", ["short", "extra", "more_points", "more_beams", "fewer_beams"])
def test_count_and_record_errors(tmp_path, change):
    path = write_data(tmp_path)
    lines = path.read_text().splitlines()
    beams = BEAMS
    if change == "short":
        lines.pop()
    elif change == "extra":
        lines.append("12 13")
    elif change == "more_points":
        lines[0] = "3 2 1"
    elif change == "more_beams":
        beams = BEAMS+((0, 0),)
    else:
        beams = BEAMS[:2]
    path.write_text("\n".join(lines))
    with pytest.raises(LegacyInputParseError, match="record"):
        parse_experimental_intensities(path, beam_indices_hk=beams)


def test_scan_shape_check_and_duplicate_scan_points(tmp_path):
    scan = LegacyScanDefinition((1, 0, 0), (0, 1, 0), (0, 0), (0, 0), (0, 0), 2, 2)
    data = parse_experimental_intensities(write_data(tmp_path), beam_indices_hk=BEAMS, scan=scan)
    assert len(generate_scan_tilts(scan)) == data.intensities.shape[0] == 4
    assert len(np.unique(generate_scan_tilts(scan), axis=0)) == 1
    for mismatch in (replace(scan, scan_columns=3), replace(scan, scan_columns=1, scan_rows=4)):
        with pytest.raises(LegacyInputParseError, match="scan dimensions"):
            parse_experimental_intensities(write_data(tmp_path), beam_indices_hk=BEAMS, scan=mismatch)


@pytest.mark.parametrize("filename", ["in.dat", "startingin.dat", "bestin.dat"])
def test_reference_files_independently(filename):
    path = REFERENCE/filename
    scan = parse_for007_scan(REFERENCE/"FOR007.DAT", experimental_header_path=path)
    data = parse_experimental_intensities(path, beam_indices_hk=REFERENCE_BEAMS, scan=scan)
    assert data.intensities.shape == data.variances.shape == data.weights.shape == (1521, 9)
    assert len(generate_scan_tilts(scan)) == 1521 and data.gradient_option == 1
    # Independent raw-row indexing, no implementation reshape used for expectations.
    lines = [line for line in path.read_text().splitlines()[1:] if line.strip()]
    assert len(lines) == 3*9*39
    for block_index, array in enumerate((data.intensities, data.variances, data.weights)):
        assert array.size == 13689 and array.dtype == np.float64
        for beam in range(9):
            for row in range(39):
                raw = [float(value) for value in lines[block_index*351+beam*39+row].split()]
                assert len(raw) == 39
                np.testing.assert_array_equal(array[row*39:(row+1)*39, beam], raw)
    assert data.intensities[0, 0] == -490.594 and data.intensities[-1, -1] == -18066.
    assert data.variances[0, 0] == 26233100 and data.variances[-1, -1] == 3608530
    assert np.count_nonzero(data.intensities < 0) == 6100
    assert np.count_nonzero(data.variances == 0) == 117
    assert np.count_nonzero(data.weights == 0) == 6140
    assert data.weights.sum() == 7549


def test_model_validation_snapshot_and_zero_beams(tmp_path):
    data = parse_experimental_intensities(write_data(tmp_path), beam_indices_hk=BEAMS)
    values = data.intensities.copy()
    copied = replace(data, intensities=values)
    values[:] = 999
    np.testing.assert_array_equal(copied.intensities, data.intensities)
    for changes in (dict(intensities=np.zeros((4, 3), dtype=np.float32)),
                    dict(weights=np.zeros((3, 4))), dict(variances=np.full((4, 3), np.nan)),
                    dict(scan_columns=True), dict(gradient_option=1.5), dict(beam_indices_hk=((True, 0),))):
        with pytest.raises((ValueError, TypeError)):
            replace(data, **changes)
    path = tmp_path/"empty-beams.dat"
    path.write_text("1 1 -5\n")
    empty = parse_experimental_intensities(path, beam_indices_hk=())
    assert empty.intensities.shape == (1, 0) and empty.gradient_option == -5


def test_wrong_scan_type_and_missing_file(tmp_path):
    with pytest.raises(TypeError, match="scan"):
        parse_experimental_intensities(tmp_path/"missing", beam_indices_hk=BEAMS, scan=1)
    with pytest.raises(LegacyInputParseError, match="readable UTF-8"):
        parse_experimental_intensities(tmp_path/"missing", beam_indices_hk=BEAMS)

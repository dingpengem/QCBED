from dataclasses import FrozenInstanceError, replace
from pathlib import Path
import math

import numpy as np
import pytest

from qcbed import LegacyScanDefinition, parse_for007_scan, generate_scan_tilts, parse_for007
from qcbed.legacy_io import LegacyInputParseError


REFERENCE = Path(__file__).parents[1]/"tests/fixtures/legacy_reference"


def definition(**changes):
    values = dict(basis_1=(1, 0, 0), basis_2=(0, 1, 0), corner_a=(-1, -2),
                  corner_b=(1, -2), corner_c=(-1, 2), scan_columns=3, scan_rows=2)
    values.update(changes)
    return LegacyScanDefinition(**values)


def write_scan(tmp_path, section=None, header=b"3 2 0\n"):
    # Reuse the protected reference prefix, with an explicitly separate header.
    lines = (REFERENCE/"FOR007.DAT").read_text().splitlines()
    stop = lines.index("in.dat   ")+1
    section = "1 0 0 0 1 0\n-1 -2\n1 -2\n-1 2\n" if section is None else section
    path = tmp_path/"scan.dat"
    path.write_text("\n".join(lines[:stop])+"\n"+section, encoding="utf-8")
    (tmp_path/"in.dat").write_bytes(header)
    return path


def test_normal_scan_column_fast_endpoints_shape_dtype():
    expected = [[-1, -2], [0, -2], [1, -2], [-1, 2], [0, 2], [1, 2]]
    tilts = generate_scan_tilts(definition())
    np.testing.assert_array_equal(tilts, expected)
    assert tilts.dtype == np.float64 and tilts.shape == (6, 2)
    assert tilts.flags.c_contiguous
    assert not np.any(np.all(tilts == 0, axis=1))  # Crossing ranges need not sample zero.


@pytest.mark.parametrize("corners,expected", [
    (((-1, 0), (1, 0), (-1, 0)), [[-1, 0], [0, 0], [1, 0]]*2),
    (((0, -1), (0, 1), (0, -1)), [[0, -1], [0, 0], [0, 1]]*2),
    (((2, 3), (4, 3), (2, 5)), [[2, 3], [3, 3], [4, 3], [2, 5], [3, 5], [4, 5]]),
    (((2, 3), (2, 3), (2, 3)), [[2, 3]]*6),
])
def test_h_only_k_only_nonzero_ranges_and_duplicate_scans(corners, expected):
    a, b, c = corners
    np.testing.assert_array_equal(generate_scan_tilts(definition(corner_a=a, corner_b=b, corner_c=c)), expected)


def test_sheared_scan_nonorthogonal_basis_and_ignored_third_components():
    scan = definition(basis_1=(2, 1, 99), basis_2=(-1, 3, -77),
        corner_a=(1, 2), corner_b=(3, 4), corner_c=(-1, 3))
    expected = [[0, 7], [1, 11], [2, 15], [-5, 8], [-4, 12], [-3, 16]]
    np.testing.assert_array_equal(generate_scan_tilts(scan), expected)
    np.testing.assert_array_equal(generate_scan_tilts(replace(scan, basis_1=(2, 1, 0), basis_2=(-1, 3, 0))), expected)


@pytest.mark.parametrize("columns,rows", [(1, 1), (1, 3), (3, 1)])
def test_singleton_preserved_but_generation_rejects_source_division_by_zero(tmp_path, columns, rows):
    path = write_scan(tmp_path, header=f"{columns} {rows} 0\n".encode())
    scan = parse_for007_scan(path)
    assert (scan.scan_columns, scan.scan_rows) == (columns, rows)
    with pytest.raises(ValueError, match="singleton.*undefined"):
        generate_scan_tilts(scan)


def test_record_continuation_and_prefix_compatibility(tmp_path):
    path = write_scan(tmp_path, "1D0,0,0\n0,1d0,0\n-1,\n-2\n1 -2\n-1 2\nREFINEMENT IS NOT PARSED /\n",
                      b"3,\n2,0\n\xff this is not a header and must not be decoded")
    parsed = parse_for007_scan(path)
    assert parsed == definition()
    assert parse_for007(path).experimental_intensity_filename == "in.dat"


def test_header_override_repetition_and_read_boundaries(tmp_path):
    path = write_scan(tmp_path, "1 2*0 0 1 0 999\n-1 -2 999\n1 -2\n-1 2\n")
    header = tmp_path/"other.dat"
    header.write_text("4 2 0\nbody is unparsed\n")
    scan = parse_for007_scan(path, experimental_header_path=header)
    assert scan.scan_columns == 4 and scan.scan_rows == 2
    assert scan.basis_1 == (1, 0, 0) and scan.corner_a == (-1, -2)


@pytest.mark.parametrize("filename", ["FOR007.DAT", "STARTINGFOR007.DAT"])
def test_real_scan_source_regression(filename):
    scan = parse_for007_scan(REFERENCE/filename)
    assert scan.basis_1 == (1, 1, 0) and scan.basis_2 == (1, -1, 0)
    assert scan.corner_a == (-.508156, .407641)
    assert scan.corner_b == (.407641, .508156)
    assert scan.corner_c == (-.407641, -.508156)
    tilts = generate_scan_tilts(scan)
    assert tilts.shape == (1521, 2) and tilts.dtype == np.float64
    # Directly evaluated h=u+v,k=u-v at A,B,C and B+C-A, plus central pixel.
    expected = [[-.100515, -.915797], [.915797, -.100515],
                [-.915797, .100515], [.100515, .915797], [0, 0]]
    np.testing.assert_allclose(tilts[[0, 38, 38*39, -1, 19*39+19]], expected, atol=3e-16, rtol=0)
    np.testing.assert_allclose(tilts[1]-tilts[0], [1.016312/38, .815282/38], atol=2e-16, rtol=0)


def test_literal_source_normalization_and_dilation_flow():
    scan = definition(basis_1=(2, 1, 0), basis_2=(1, -3, 0),
                      corner_a=(-.3, .2), corner_b=(.7, .4), corner_c=(-.2, .9))
    # Independent scalar execution of ONEIT's steps with a nontrivial norm
    # computed by the driver's sign branches at reciprocal gamma=70 degrees.
    lx = math.sqrt(5-4*math.cos(3.1415926535*110/180))
    ly = math.sqrt(10-6*math.cos(3.1415926535*70/180))
    norm, dilation = ly/lx, 1.17
    a, b, c = [[p[0]/dilation, p[1]*norm/dilation] for p in
               (scan.corner_a, scan.corner_b, scan.corner_c)]
    expected = []
    for row in range(2):
        for col in range(3):
            u = a[0]+row*(c[0]-a[0])+col*(b[0]-a[0])/2
            v = (a[1]+col*(b[1]-a[1])/2+row*(c[1]-a[1]))/norm
            expected.append([2*u+v, u-3*v])
    np.testing.assert_allclose(generate_scan_tilts(scan)/dilation, expected, atol=9e-16, rtol=0)


@pytest.mark.parametrize("name,bad", [("scan_columns", 0), ("scan_rows", -1),
    ("scan_columns", True), ("scan_rows", 2.), ("basis_1", (0, 0, 1)),
    ("basis_2", (0, 0, 0)), ("basis_1", (1, 0)), ("corner_a", (np.nan, 0)),
    ("corner_b", (0, np.inf)), ("corner_c", (True, 0)), ("corner_a", "12")])
def test_model_validation(name, bad):
    with pytest.raises((ValueError, TypeError), match=name):
        definition(**{name: bad})


@pytest.mark.parametrize("header", [b"0 2 0", b"1.0 2 0", b"2 NaN 0", b"2 2 word", b"2"])
def test_header_errors_have_context(tmp_path, header):
    path = write_scan(tmp_path, header=header)
    with pytest.raises(LegacyInputParseError, match="in.dat.*scan dimensions header.*record"):
        parse_for007_scan(path)


@pytest.mark.parametrize("section,context", [("0 0 1\n0 1 0\n", "scan basis.*record 36"),
    ("NaN 0 0 0 1 0\n", "scan basis"), ("1 0 0 0 1 0\n1+ 0\n", "scan corner_a"),
    ("1 0 0 0 1 0\n1 2\n3 /\n", "scan corner_b"),
    ("1 0 0 0 1 0\n1 2\n3 4\n", "scan corner_c")])
def test_scan_errors_have_context(tmp_path, section, context):
    path = write_scan(tmp_path, section)
    with pytest.raises(LegacyInputParseError, match=f"scan.dat.*{context}") as error:
        parse_for007_scan(path)
    assert "record" in str(error.value) and "expected" in str(error.value) and "observed" in str(error.value)


def test_missing_header_no_inferred_counts(tmp_path):
    path = write_scan(tmp_path)
    with pytest.raises(LegacyInputParseError, match="scan dimensions header.*missing.dat"):
        parse_for007_scan(path, experimental_header_path=tmp_path/"missing.dat")


def test_immutability_and_generated_overflow():
    a = [1, 2]
    scan = definition(corner_a=a)
    a[0] = 999
    assert scan.corner_a == (1, 2)
    with pytest.raises(FrozenInstanceError):
        scan.scan_columns = 10
    with pytest.raises(ValueError, match="finite"):
        generate_scan_tilts(definition(corner_a=(-1e308, 0), corner_b=(1e308, 0)))
    with pytest.raises(TypeError, match="LegacyScanDefinition"):
        generate_scan_tilts(None)

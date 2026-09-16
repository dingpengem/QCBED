from pathlib import Path

import numpy as np
import pytest

from qcbed.legacy_io import LegacyInputParseError, parse_for007


REFERENCE_FOR007 = Path(__file__).parents[1] / "reference_data" / "FOR007.DAT"


def test_parse_small_reference_for007() -> None:
    parsed = parse_for007(REFERENCE_FOR007)

    # Source: oneitpara.f READ(7,*) KV,KVOPT,INDZ; FOR0013.DAT logs 200.00000 KV
    # and zone axis 0 0 1 for this reference run.
    assert parsed.title == "PACBED EMPAD MONASH E.T. Trondheim"
    assert parsed.energy_kv == 200.0
    assert parsed.energy_refined is False
    assert parsed.zone_axis == (0, 0, 1)

    # Source: oneitpara.f READ(7,*) NSLICE,MESHX,MESHY; the log reports 1 slice
    # and 16 by 16 sampling points in X and Y.
    assert parsed.number_of_slices == 1
    assert parsed.mesh.mesh_x == 16
    assert parsed.mesh.mesh_y == 16

    # Source: oneitpara.f READ(7,*) CELL2D(...),DELTAZ,XSHIFT,YSHIFT; the log
    # repeats this first and only slice geometry record.
    first_geometry = parsed.slices[0]
    last_geometry = parsed.slices[-1]
    assert first_geometry == last_geometry
    assert first_geometry.cell_a_angstrom == pytest.approx(2.8630754)
    assert first_geometry.cell_b_angstrom == pytest.approx(2.8630754)
    assert first_geometry.lattice_gamma_degrees == pytest.approx(90.0)
    assert first_geometry.thickness_angstrom == pytest.approx(2.0245)
    assert first_geometry.x_shift_angstrom == pytest.approx(0.0)
    assert first_geometry.y_shift_angstrom == pytest.approx(0.0)
    assert np.asarray(parsed.slices, dtype=object).shape == (1,)
    assert np.asarray(parsed.atoms_by_slice, dtype=object).shape == (1, 2)

    # Source: oneitpara.f reads LABEL then ZT for each species; SFMAP reads
    # NBALL/RTYPE and the two S atom records from the following records.
    assert parsed.species[0].label == "AL"
    assert parsed.species[0].atomic_number == 13
    assert len(parsed.atoms_by_slice) == 1
    assert len(parsed.atoms_by_slice[0]) == 2
    assert parsed.atoms_by_slice[0][0].occupancy == pytest.approx(0.5)
    assert parsed.atoms_by_slice[0][1].x == pytest.approx(0.5)

    # Source: READ(7,452) BMNAME and READ(7,452) EXPINT; records are out.dat
    # and in.dat in the actual reference input.
    assert parsed.output_filename == "out.dat"
    assert parsed.experimental_intensity_filename == "in.dat"


def test_parser_preserves_float64_compatible_values() -> None:
    parsed = parse_for007(REFERENCE_FOR007)
    geometry = parsed.slices[0]
    values = np.asarray(
        [
            parsed.energy_kv,
            geometry.cell_a_angstrom,
            geometry.cell_b_angstrom,
            geometry.lattice_gamma_degrees,
            geometry.thickness_angstrom,
        ],
        dtype=np.float64,
    )
    assert values.dtype == np.float64
    assert np.all(np.isfinite(values))


def test_truncated_input_reports_section_and_record(tmp_path: Path) -> None:
    truncated = tmp_path / "truncated.FOR007.DAT"
    truncated.write_text("title\n200 0 0 0 1\n", encoding="utf-8")

    with pytest.raises(LegacyInputParseError, match=r"truncated\.FOR007\.DAT.*slice count"):
        parse_for007(truncated)


def test_invalid_integer_reports_observed_content(tmp_path: Path) -> None:
    lines = REFERENCE_FOR007.read_text(encoding="utf-8").splitlines()
    lines[2] = "1 sixteen 16"
    invalid = tmp_path / "invalid.FOR007.DAT"
    invalid.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(LegacyInputParseError, match=r"invalid\.FOR007\.DAT.*mesh_x.*sixteen"):
        parse_for007(invalid)


def test_invalid_slice_count_reports_expected_positive_count(tmp_path: Path) -> None:
    lines = REFERENCE_FOR007.read_text(encoding="utf-8").splitlines()
    lines[2] = "0 16 16"
    invalid = tmp_path / "invalid-slices.FOR007.DAT"
    invalid.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(LegacyInputParseError, match=r"invalid-slices\.FOR007\.DAT.*slice count"):
        parse_for007(invalid)
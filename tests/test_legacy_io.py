from pathlib import Path

import numpy as np
import pytest

from qcbed.legacy_io import LegacyInputParseError, parse_for007
from qcbed import AtomSite, SliceGeometry


REFERENCE_FOR007 = Path(__file__).parents[1] / "tests/fixtures/legacy_reference" / "FOR007.DAT"


@pytest.fixture
def synthetic_input() -> list[str]:
    # Records follow PROGRAM MULTIS/SFMAP reads, with distinct a/b and x/y.
    return """Synthetic fractional coordinates
200 0 0 0 1
1 6 4
3 5 80 2 0.25 -0.375
2
AL
13
SI
14
2.5
0.5
1
10 0
1
0
0.2 0
0 0
NO
0.8 0.9
0 0
0 0 0 1
1 S
2 0.125 1 0.75 0 0.6
out.dat
1
0 0
in.dat
""".splitlines()


@pytest.mark.parametrize("record_type", ["S", "G"])
def test_synthetic_fractional_shifts_and_atom_semantics(
    tmp_path: Path, synthetic_input: list[str], record_type: str
) -> None:
    if record_type == "G":
        synthetic_input[21:23] = ["1 G", "2 0.125 0.75 0.6"]
    path = tmp_path / "synthetic.dat"
    path.write_text("\n".join(synthetic_input), encoding="utf-8")
    parsed = parse_for007(path)
    geometry = parsed.slices[0]
    assert isinstance(geometry, SliceGeometry)
    assert geometry.x_shift_fractional == 0.25
    assert geometry.y_shift_fractional == -0.375
    assert geometry.cell_a_angstrom == 3.0
    assert geometry.cell_b_angstrom == 5.0
    atom = parsed.atoms_by_slice[0][0]
    assert isinstance(atom, AtomSite)
    assert atom.species_index_one_based == 2
    assert atom.x_fractional == 0.125
    assert atom.y_fractional == 0.75
    assert atom.starting_occupancy == 0.6  # Not multiplied by (1 - vacancy).
    assert atom.record_type == record_type
    assert atom.x_mask == (1 if record_type == "S" else 0)
    assert atom.y_mask == 0
    assert atom.species_index_zero_based(len(parsed.species)) == 1
    assert parsed.species[atom.species_index_zero_based(2)].label == "SI"


@pytest.mark.parametrize("index", [-1, 0, 3])
@pytest.mark.parametrize("record_type", ["S", "G"])
def test_parser_rejects_out_of_range_species(
    tmp_path: Path, synthetic_input: list[str], index: int, record_type: str
) -> None:
    synthetic_input[21] = f"1 {record_type}"
    synthetic_input[22] = (
        f"{index} 0.125 1 0.75 0 0.6" if record_type == "S"
        else f"{index} 0.125 0.75 0.6"
    )
    path = tmp_path / "invalid-species.dat"
    path.write_text("\n".join(synthetic_input), encoding="utf-8")
    with pytest.raises(LegacyInputParseError, match=r"atom 1.*record 23.*1\.\.2"):
        parse_for007(path)


@pytest.mark.parametrize("index, expected", [(1, 0), (2, 1)])
def test_species_conversion_checks_endpoints(index: int, expected: int) -> None:
    assert AtomSite(index, 0.0, 0.0, 1.0, "S").species_index_zero_based(2) == expected


@pytest.mark.parametrize("index", [-1, 0, 3])
def test_species_conversion_rejects_invalid_index(index: int) -> None:
    with pytest.raises(ValueError, match="species_index_one_based"):
        AtomSite(index, 0.0, 0.0, 1.0, "S").species_index_zero_based(2)


@pytest.mark.parametrize("count", [0, -1])
def test_species_conversion_rejects_invalid_count(count: int) -> None:
    with pytest.raises(ValueError, match="species_count"):
        AtomSite(1, 0.0, 0.0, 1.0, "S").species_index_zero_based(count)


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
    assert first_geometry.x_shift_fractional == pytest.approx(0.0)
    assert first_geometry.y_shift_fractional == pytest.approx(0.0)
    assert np.asarray(parsed.slices, dtype=object).shape == (1,)
    assert np.asarray(parsed.atoms_by_slice, dtype=object).shape == (1, 2)

    # Source: oneitpara.f reads LABEL then ZT for each species; SFMAP reads
    # NBALL/RTYPE and the two S atom records from the following records.
    assert parsed.species[0].label == "AL"
    assert parsed.species[0].atomic_number == 13
    assert len(parsed.atoms_by_slice) == 1
    assert len(parsed.atoms_by_slice[0]) == 2
    assert parsed.atoms_by_slice[0][0].starting_occupancy == pytest.approx(0.5)
    assert parsed.atoms_by_slice[0][1].x_fractional == pytest.approx(0.5)

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


@pytest.fixture
def parse_lines(tmp_path: Path):
    def parse(lines):
        path = tmp_path / "syntax.FOR007.DAT"
        path.write_text("\n".join(lines), encoding="utf-8")
        return parse_for007(path)
    return parse


@pytest.mark.parametrize("bad", ["NaN", "Inf", "+Infinity", "-inf", "1e999"])
@pytest.mark.parametrize("record, field", [
    (1, 0), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5),
    (9, 0), (10, 0), (12, 0), (15, 0), (16, 0), (18, 0),
    (20, 1), (22, 1), (22, 3), (22, 5),
])
def test_nonfinite_scientific_fields_rejected(synthetic_input, parse_lines, bad, record, field):
    tokens = synthetic_input[record].split()
    tokens[field] = bad
    synthetic_input[record] = ' '.join(tokens)
    with pytest.raises(LegacyInputParseError) as exc:
        parse_lines(synthetic_input)
    message = str(exc.value)
    assert 'syntax.FOR007.DAT' in message
    assert "section '" in message
    assert f'record {record + 1}:' in message
    assert bad in message
    assert 'finite' in message


@pytest.mark.parametrize("bad", ["0.5", "1.0", "NaN", "text"])
@pytest.mark.parametrize("record, field", [
    (1, 1), (1, 2), (2, 0), (2, 1), (4, 0), (6, 0),
    (11, 0), (12, 1), (13, 0), (14, 0), (15, 1), (16, 1),
    (19, 0), (20, 0), (21, 0), (22, 0), (22, 2), (22, 4),
    (24, 0), (25, 0),
])
def test_integer_fields_never_coerce_float(synthetic_input, parse_lines, bad, record, field):
    tokens = synthetic_input[record].split()
    tokens[field] = bad
    synthetic_input[record] = ' '.join(tokens)
    with pytest.raises(LegacyInputParseError) as exc:
        parse_lines(synthetic_input)
    assert 'integer' in str(exc.value)
    assert f'record {record + 1}:' in str(exc.value)
    assert bad in str(exc.value)


@pytest.mark.parametrize("record, field, bad", [
    (1, 1, '2'), (12, 1, '-1'), (14, 0, '3'), (15, 1, '2'),
    (16, 1, '-1'), (19, 1, '2'), (20, 0, '3'), (22, 2, '2'),
    (22, 4, '-1'),
])
def test_mask_and_flag_domains(synthetic_input, parse_lines, record, field, bad):
    tokens = synthetic_input[record].split()
    tokens[field] = bad
    synthetic_input[record] = ' '.join(tokens)
    with pytest.raises(LegacyInputParseError, match='integer in'):
        parse_lines(synthetic_input)


@pytest.mark.parametrize("record, field", [(1, 0), (3, 0), (3, 1), (3, 3), (12, 0)])
@pytest.mark.parametrize("bad", ['0', '-1'])
def test_positive_scientific_values(synthetic_input, parse_lines, record, field, bad):
    tokens = synthetic_input[record].split()
    tokens[field] = bad
    synthetic_input[record] = ' '.join(tokens)
    with pytest.raises(LegacyInputParseError, match='positive'):
        parse_lines(synthetic_input)


@pytest.mark.parametrize('gamma', ['0', '180', '-90'])
def test_model_failure_keeps_parser_context(synthetic_input, parse_lines, gamma):
    synthetic_input[3] = f'3 5 {gamma} 2 .25 -.375'
    with pytest.raises(LegacyInputParseError) as exc:
        parse_lines(synthetic_input)
    message = str(exc.value)
    assert "syntax.FOR007.DAT: section 'slice geometry 1', record 4:" in message
    assert synthetic_input[3] in message
    assert 'between 0 and 180' in message


def test_numeric_forms_and_repetition(synthetic_input, parse_lines):
    synthetic_input[1] = '2D2, 3*0, +1'
    synthetic_input[3] = '3., 5d0, 8e1, +2, +.25, -.375'
    synthetic_input[18] = '2*8.5d-1'
    synthetic_input[19] = '2*0'
    synthetic_input[22] = '2,1.25-1,1,7.5E-1,0,6d-1'
    parsed = parse_lines(synthetic_input)
    assert parsed.energy_kv == 200
    assert parsed.zone_axis == (0, 0, 1)
    assert parsed.slices[0].x_shift_fractional == .25
    assert parsed.slices[0].y_shift_fractional == -.375
    atom = parsed.atoms_by_slice[0][0]
    assert (atom.x_fractional, atom.y_fractional, atom.starting_occupancy) == (.125, .75, .6)


@pytest.mark.parametrize('literal', ['1+', '1-', '0.5+', '2.0D3-'])
@pytest.mark.parametrize('record, field', [(1, 2), (3, 4)])
def test_trailing_signs_rejected_in_integer_and_real_fields(
    synthetic_input, parse_lines, literal, record, field
):
    tokens = synthetic_input[record].split()
    tokens[field] = literal
    synthetic_input[record] = ' '.join(tokens)
    with pytest.raises(LegacyInputParseError) as exc:
        parse_lines(synthetic_input)
    assert f'record {record + 1}:' in str(exc.value)
    assert repr(literal) in str(exc.value)


@pytest.mark.parametrize('literal, expected', [
    ('+1', 1.0), ('-1', -1.0), ('+0.5', 0.5),
    ('1E+3', 1000.0), ('1e-3', 0.001),
    ('2.0D+4', 20000.0), ('2.0d-4', 0.0002),
])
def test_leading_and_exponent_signs_in_real_fields(
    synthetic_input, parse_lines, literal, expected
):
    synthetic_input[3] = f'3 5 80 2 {literal} -.375'
    assert parse_lines(synthetic_input).slices[0].x_shift_fractional == expected


@pytest.mark.parametrize('literal, expected', [('+1', 1), ('-1', -1)])
def test_leading_signs_in_integer_fields(synthetic_input, parse_lines, literal, expected):
    synthetic_input[1] = f'200 0 {literal} 0 1'
    assert parse_lines(synthetic_input).zone_axis == (expected, 0, 1)


@pytest.mark.parametrize('beam_records', [
    ['2, 0, 0, -1, +1'], ['2, 0', '0, -1', '+1'], ['2', '0 0', '-1 1'],
])
def test_beam_count_and_indices_share_one_read(synthetic_input, parse_lines, beam_records):
    synthetic_input[24:26] = beam_records
    assert parse_lines(synthetic_input).experimental_intensity_filename == 'in.dat'


def test_reads_continue_without_flattening_character_records(synthetic_input, parse_lines):
    synthetic_input[23] = ' output name with spaces.dat  '
    synthetic_input[22:23] = ['2 .125', '1 .75', '0 .6']
    synthetic_input[3:4] = ['3, 5,', '', '80 2', '.25 -.375']
    synthetic_input[1:2] = ['200', '0 0', '0 1']
    parsed = parse_lines(synthetic_input)
    assert parsed.slices[0].y_shift_fractional == -.375
    assert parsed.output_filename_raw == ' output name with spaces.dat  '
    assert parsed.output_filename == ' output name with spaces.dat'
    assert parsed.atoms_by_slice[0][0].starting_occupancy == .6


@pytest.mark.parametrize('raw', ["'S'", '\"S\"', "1*'S'", 'S'])
def test_quoted_list_directed_record_types(synthetic_input, parse_lines, raw):
    synthetic_input[21] = f'1, {raw}'
    atom = parse_lines(synthetic_input).atoms_by_slice[0][0]
    assert atom.record_type == 'S'
    assert atom.record_type_raw == raw


@pytest.mark.parametrize('raw', ['s', "'s'", "' S'", "'S '", 'SOLID'])
def test_record_type_not_normalized(synthetic_input, parse_lines, raw):
    synthetic_input[21] = f'1 {raw}'
    with pytest.raises(LegacyInputParseError, match='exactly S or G'):
        parse_lines(synthetic_input)


@pytest.mark.parametrize('raw, expected', [
    ('Al', 'Al'), ('C ', 'C '), (' C', ' C'), ('C', 'C '),
    ('ALUMINUM', 'AL'), ("'AL'", "'A"), ('\"AL\"', '\"A'),
])
def test_species_labels_use_formatted_a2_not_list_quoting(synthetic_input, parse_lines, raw, expected):
    synthetic_input[5] = raw
    parsed = parse_lines(synthetic_input)
    assert parsed.species[0].label == expected
    assert parsed.species[0].label_raw == raw


@pytest.mark.parametrize('raw, enabled', [
    ('YES_ABSORPTION  ', True), ('YEanything', True), ('yes', False),
    (' YES', False), ("'YES'", False), ('NO', False), ('', False),
])
def test_absorption_uses_literal_first_two_characters(synthetic_input, parse_lines, raw, enabled):
    synthetic_input[17] = raw
    if enabled:
        synthetic_input[18:18] = ['0 0 1 1 1 -1', '0 0 1 0 1 1']
    parsed = parse_lines(synthetic_input)
    assert parsed.absorption_raw == raw
    assert parsed.absorption_enabled is enabled


@pytest.mark.parametrize('field, bad', [(0, 'nan'), (1, '-1'), (3, '2'), (5, '-2'), (5, '1.0')])
def test_absorption_coefficients_and_masks_validated(synthetic_input, parse_lines, field, bad):
    synthetic_input[17] = 'YES'
    values = '0 0 1 1 1 -1'.split()
    values[field] = bad
    synthetic_input[18:18] = [' '.join(values), '0 0 1 0 1 0']
    with pytest.raises(LegacyInputParseError, match='absorption coefficients 1.*record 19'):
        parse_lines(synthetic_input)


@pytest.mark.parametrize('bad', [
    '200,,0,0,0,1', ',200 0 0 0 1', '200 /', '200 2* 0 1',
    '0*200 0 0 0 1', '200 ! comment', '(200,0) 0 0 0 1',
    '200 0 0 0 1 /',
])
def test_unsupported_syntax_is_explicit(synthetic_input, parse_lines, bad):
    synthetic_input[1] = bad
    with pytest.raises(LegacyInputParseError) as exc:
        parse_lines(synthetic_input)
    assert 'record 2:' in str(exc.value)
    assert 'expected' in str(exc.value)
    assert 'observed' in str(exc.value)


def test_unterminated_quote(synthetic_input, parse_lines):
    synthetic_input[21] = "1 'S"
    with pytest.raises(LegacyInputParseError, match='closed quoted character'):
        parse_lines(synthetic_input)


def test_continuation_error_reports_actual_offending_record(synthetic_input, parse_lines):
    synthetic_input[1:2] = ['200 0', '0 bad 1']
    with pytest.raises(LegacyInputParseError, match="zone axis.*record 3:.*integer.*bad"):
        parse_lines(synthetic_input)


def test_extra_values_are_discarded_at_end_of_read(synthetic_input, parse_lines):
    synthetic_input[4] = '2 999'
    assert len(parse_lines(synthetic_input).species) == 2


def test_multiple_regions_slices_and_species(parse_lines):
    # Explicit READ-by-READ fixture: two distinct geometries and S/G records.
    lines = """Two regions
200 0 0 0 1
2 6 4
3 5 80 2 .25 -.375
4 6 90 3 -.5 .125
2
Al
13
Si
14
2.5
.5
2
10 0
12 1
1 1
1
.1 0 .2 1
2
.3 1 .4 0
0 0
NO
.8 .9
0 1
1 .1 .2 .3
1.1 1.2
1 0
2 .4 .5 .6
1 'S'
1 .1 0 .2 1 .7
1 "G"
2 .3 .4 .8
out.dat
2 0 0 -1 1
in.dat
""".splitlines()
    parsed = parse_lines(lines)
    assert parsed.number_of_slices == 2
    assert [s.label for s in parsed.species] == ['Al', 'Si']
    assert [s.thickness_angstrom for s in parsed.slices] == [2, 3]
    first, second = [atoms[0] for atoms in parsed.atoms_by_slice]
    assert (first.species_index_one_based, first.record_type, first.starting_occupancy) == (1, 'S', .7)
    assert (second.species_index_one_based, second.record_type, second.starting_occupancy) == (2, 'G', .8)
    assert (second.x_mask, second.y_mask) == (0, 0)
    assert parsed.experimental_intensity_filename == 'in.dat'


def test_region_slice_counts_must_cover_geometry(synthetic_input, parse_lines):
    synthetic_input[13] = '2'
    with pytest.raises(LegacyInputParseError, match='region slice counts.*summing to slice count 1'):
        parse_lines(synthetic_input)


def test_geometry_model_error_on_continued_read(synthetic_input, parse_lines):
    synthetic_input[3:4] = ['3 5 180', '2 .25 -.375']
    with pytest.raises(LegacyInputParseError, match="slice geometry 1.*record 4:.*180.*3 5 180"):
        parse_lines(synthetic_input)


def test_raw_title_and_filename_records(synthetic_input, parse_lines):
    synthetic_input[0] = '  title  '
    synthetic_input[26] = '  a long experimental filename.dat  '
    parsed = parse_lines(synthetic_input)
    assert parsed.title_raw == '  title  '
    assert parsed.experimental_intensity_filename_raw == '  a long experimental filename.dat  '
    assert parsed.experimental_intensity_filename == '  a long experimental filename.dat'

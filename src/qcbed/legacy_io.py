"""Parser for the confirmed first-stage FOR007.DAT records."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Callable, Sequence, TypeVar

from .models import MeshConfig, SliceGeometry


T = TypeVar("T")
_INT_RE = re.compile(r"^[+-]?\d+$")


class LegacyInputParseError(ValueError):
    """A source-record-aware FOR007 parsing error."""


@dataclass(frozen=True)
class Species:
    """One source-confirmed atomic species entry."""

    label: str
    atomic_number: int


@dataclass(frozen=True)
class AtomSite:
    """One source-confirmed atom record from ``SFMAP``."""

    species_index: int
    x: float
    y: float
    occupancy: float
    x_mask: int = 0
    y_mask: int = 0


@dataclass(frozen=True)
class LegacyInput:
    """The confirmed first-stage subset of a legacy FOR007 input."""

    title: str
    energy_kv: float
    energy_refined: bool
    zone_axis: tuple[int, int, int]
    mesh: MeshConfig
    slices: tuple[SliceGeometry, ...]
    species: tuple[Species, ...]
    atoms_by_slice: tuple[tuple[AtomSite, ...], ...]
    output_filename: str
    experimental_intensity_filename: str

    @property
    def number_of_slices(self) -> int:
        return len(self.slices)


class _Records:
    def __init__(self, filename: Path) -> None:
        self.filename = filename
        self.lines = filename.read_text(encoding="utf-8").splitlines()
        self.position = 0

    @property
    def record_number(self) -> int:
        return self.position + 1

    def next(self, section: str) -> str:
        if self.position >= len(self.lines):
            raise LegacyInputParseError(
                f"{self.filename}: section '{section}', record {self.record_number}: "
                "unexpected end of file"
            )
        line = self.lines[self.position]
        self.position += 1
        return line

    def error(self, section: str, observed: str, expected: str) -> LegacyInputParseError:
        return LegacyInputParseError(
            f"{self.filename}: section '{section}', record {self.position}: "
            f"expected {expected}; observed {observed!r}"
        )


def _parse_int(records: _Records, section: str, token: str) -> int:
    if not _INT_RE.fullmatch(token):
        raise records.error(section, token, "an integer")
    return int(token)


def _parse_float(records: _Records, section: str, token: str) -> float:
    try:
        return float(token)
    except ValueError as exc:
        raise records.error(section, token, "a floating-point value") from exc


def _tokens(records: _Records, section: str, expected: int) -> list[str]:
    line = records.next(section)
    tokens = line.split()
    if len(tokens) != expected:
        raise records.error(section, line, f"{expected} values")
    return tokens


def _values(
    records: _Records,
    section: str,
    count: int,
    converter: Callable[[_Records, str, str], T],
) -> list[T]:
    values: list[T] = []
    while len(values) < count:
        line = records.next(section)
        tokens = line.split()
        if not tokens:
            raise records.error(section, line, f"{count} values in list-directed input")
        for token in tokens:
            values.append(converter(records, section, token))
            if len(values) == count:
                break
    return values


def _text_record(records: _Records, section: str) -> str:
    value = records.next(section).strip()
    if not value:
        raise records.error(section, value, "a non-empty character record")
    return value


def parse_for007(path: str | Path) -> LegacyInput:
    """Parse the confirmed first-stage records from a formatted FOR007 file.

    The parser follows the record order in ``PROGRAM MULTIS`` and ``SFMAP``.
    It consumes intervening confirmed records only to reach the atom and
    filename records; it does not expose or calculate their physics values.
    Later refinement and scan sections remain intentionally unsupported.
    """
    filename = Path(path)
    records = _Records(filename)

    title = records.next("title").rstrip()
    header = _tokens(records, "beam voltage and zone axis", 5)
    energy_kv = _parse_float(records, "beam voltage", header[0])
    kvopt = _parse_int(records, "energy refinement flag", header[1])
    zone_axis = tuple(_parse_int(records, "zone-axis indices", token) for token in header[2:])
    if len(zone_axis) != 3:
        raise records.error("zone-axis indices", " ".join(header[2:]), "three integers")

    mesh_tokens = _tokens(records, "slice count and mesh dimensions", 3)
    number_of_slices = _parse_int(records, "slice count", mesh_tokens[0])
    if number_of_slices <= 0:
        raise records.error("slice count", mesh_tokens[0], "a positive integer")
    mesh_x = _parse_int(records, "mesh_x", mesh_tokens[1])
    mesh_y = _parse_int(records, "mesh_y", mesh_tokens[2])
    mesh = MeshConfig(mesh_x, mesh_y)

    slices: list[SliceGeometry] = []
    for index in range(number_of_slices):
        section = f"slice geometry {index + 1}"
        values = _tokens(records, section, 6)
        slices.append(
            SliceGeometry(
                cell_a_angstrom=_parse_float(records, section, values[0]),
                cell_b_angstrom=_parse_float(records, section, values[1]),
                lattice_gamma_degrees=_parse_float(records, section, values[2]),
                thickness_angstrom=_parse_float(records, section, values[3]),
                x_shift_angstrom=_parse_float(records, section, values[4]),
                y_shift_angstrom=_parse_float(records, section, values[5]),
            )
        )

    species_count = _parse_int(records, "species count", _tokens(records, "species count", 1)[0])
    if species_count <= 0:
        raise records.error("species count", str(species_count), "a positive integer")
    species: list[Species] = []
    for index in range(species_count):
        label = _text_record(records, f"species {index + 1} label")
        atomic_number = _parse_int(
            records,
            f"species {index + 1} atomic number",
            _tokens(records, f"species {index + 1} atomic number", 1)[0],
        )
        species.append(Species(label, atomic_number))

    # These records are consumed exactly as the confirmed list-directed reads
    # in oneitpara.f; their physics interpretation is intentionally deferred.
    _parse_float(records, "DSAP", _tokens(records, "DSAP", 1)[0])
    _parse_float(records, "DSMAX", _tokens(records, "DSMAX", 1)[0])
    region_count = _parse_int(records, "region count", _tokens(records, "region count", 1)[0])
    if region_count <= 0:
        raise records.error("region count", str(region_count), "a positive integer")
    _values(records, "region depths and masks", 2 * region_count, _parse_float_or_int)
    _values(records, "region slice counts", region_count, _parse_int_token)

    for region in range(region_count):
        concflag = _parse_int_token(records, f"region {region + 1} vacancy flag", _tokens(records, f"region {region + 1} vacancy flag", 1)[0])
        expected = 4 if concflag > 0 else 2
        _tokens(records, f"region {region + 1} vacancy parameters", expected)

    _tokens(records, "lattice dilation", 2)
    absorption = _text_record(records, "absorption flag")
    if absorption[:2].upper() == "YE":
        for species_index in range(species_count):
            _tokens(records, f"absorption coefficients {species_index + 1}", 6)
    for region in range(region_count):
        _values(records, f"Debye-Waller factors region {region + 1}", species_count, _parse_float_token)
        _values(records, f"Debye-Waller masks region {region + 1}", species_count, _parse_int_token)
        _tokens(records, f"free-surface parameters region {region + 1}", 4)

    atoms_by_slice: list[tuple[AtomSite, ...]] = []
    for slice_index in range(number_of_slices):
        section = f"slice {slice_index + 1} atom header"
        header_tokens = _tokens(records, section, 2)
        atom_count = _parse_int(records, section, header_tokens[0])
        if atom_count < 0:
            raise records.error(section, header_tokens[0], "a non-negative integer atom count")
        record_type = header_tokens[1].upper()
        if record_type not in {"S", "G"}:
            raise records.error(section, header_tokens[1], "record type S or G")
        atoms: list[AtomSite] = []
        for atom_index in range(atom_count):
            atom_section = f"slice {slice_index + 1} atom {atom_index + 1}"
            if record_type == "S":
                values = _tokens(records, atom_section, 6)
                atoms.append(
                    AtomSite(
                        species_index=_parse_int(records, atom_section, values[0]),
                        x=_parse_float(records, atom_section, values[1]),
                        x_mask=_parse_int(records, atom_section, values[2]),
                        y=_parse_float(records, atom_section, values[3]),
                        y_mask=_parse_int(records, atom_section, values[4]),
                        occupancy=_parse_float(records, atom_section, values[5]),
                    )
                )
            else:
                values = _tokens(records, atom_section, 4)
                atoms.append(
                    AtomSite(
                        species_index=_parse_int(records, atom_section, values[0]),
                        x=_parse_float(records, atom_section, values[1]),
                        y=_parse_float(records, atom_section, values[2]),
                        occupancy=_parse_float(records, atom_section, values[3]),
                    )
                )
        atoms_by_slice.append(tuple(atoms))

    output_filename = _text_record(records, "output filename BMNAME")
    beam_header = _tokens(records, "output beam count", 1)
    beam_count = _parse_int(records, "output beam count", beam_header[0])
    if beam_count < 0:
        raise records.error("output beam count", beam_header[0], "a non-negative integer")
    _values(records, "output beam indices", 2 * beam_count, _parse_int_token)
    experimental_filename = _text_record(records, "experimental intensity filename EXPINT")

    return LegacyInput(
        title=title,
        energy_kv=energy_kv,
        energy_refined=kvopt == 1,
        zone_axis=zone_axis,
        mesh=mesh,
        slices=tuple(slices),
        species=tuple(species),
        atoms_by_slice=tuple(atoms_by_slice),
        output_filename=output_filename,
        experimental_intensity_filename=experimental_filename,
    )


def _parse_int_token(records: _Records, section: str, token: str) -> int:
    return _parse_int(records, section, token)


def _parse_float_token(records: _Records, section: str, token: str) -> float:
    return _parse_float(records, section, token)


def _parse_float_or_int(records: _Records, section: str, token: str) -> float | int:
    try:
        return _parse_int(records, section, token)
    except LegacyInputParseError:
        return _parse_float(records, section, token)

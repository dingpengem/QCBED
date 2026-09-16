"""Parser for the confirmed first-stage FOR007.DAT records."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import math
from collections.abc import Iterator
from typing import Callable, TypeVar

from .models import MeshConfig, SliceGeometry


T = TypeVar("T")
_INT_RE = re.compile(r"^[+-]?[0-9]+$")


class LegacyInputParseError(ValueError):
    """A source-record-aware FOR007 parsing error."""


@dataclass(frozen=True)
class Species:
    """One source-confirmed atomic species entry."""

    label: str
    atomic_number: int
    label_raw: str = ""


@dataclass(frozen=True)
class AtomSite:
    """Input S/G record: fractional positions and unadjusted STARTOCC.

    Gas randomization and vacancy corrections are not applied by the parser.
    """

    species_index_one_based: int
    x_fractional: float
    y_fractional: float
    starting_occupancy: float
    record_type: str
    x_mask: int = 0
    y_mask: int = 0
    record_type_raw: str = ""

    def species_index_zero_based(self, species_count: int) -> int:
        """Return a Python index only after checking both species bounds."""
        if isinstance(species_count, bool) or not isinstance(species_count, int):
            raise TypeError("species_count must be an integer")
        if species_count <= 0:
            raise ValueError("species_count must be positive")
        index = self.species_index_one_based
        if isinstance(index, bool) or not isinstance(index, int):
            raise TypeError("species_index_one_based must be an integer")
        if not 1 <= index <= species_count:
            raise ValueError(f"species_index_one_based must be in 1..{species_count}")
        return index - 1


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
    title_raw: str = ""
    absorption_raw: str = ""
    output_filename_raw: str = ""
    experimental_intensity_filename_raw: str = ""

    @property
    def absorption_enabled(self) -> bool:
        return self.absorption_raw[:2] == "YE"

    @property
    def number_of_slices(self) -> int:
        return len(self.slices)


@dataclass(frozen=True)
class _Token:
    raw: str
    value: str
    record: int
    quoted: bool = False
    repeat: int = 1


class _Records:
    def __init__(self, filename: Path) -> None:
        self.filename = filename
        self.position = 0
        try:
            self.lines = filename.read_text(encoding="utf-8").splitlines()
        except (OSError, UnicodeError) as exc:
            raise self.error("input", str(filename), "a readable UTF-8 file", 1) from exc

    def error(self, section: str, observed: str, expected: str,
              record: int | None = None) -> LegacyInputParseError:
        return LegacyInputParseError(
            f"{self.filename}: section '{section}', record "
            f"{self.position if record is None else record}: "
            f"expected {expected}; observed {observed!r}"
        )

    def next(self, section: str) -> str:
        if self.position == len(self.lines):
            raise self.error(section, "<EOF>", "another input record",
                             self.position + 1)
        line = self.lines[self.position]
        self.position += 1
        return line

    def read(self, section: str) -> _Read:
        # One cursor per actual Fortran READ, never one global token stream.
        return _Read(self, section)


def _lex(records: _Records, section: str, line: str) -> Iterator[_Token]:
    """Tokenize one record; preserve quoting and expand repetitions lazily."""
    i = 0
    have_value = False
    comma_pending = False
    while i < len(line):
        if line[i].isspace():
            i += 1
            continue
        if line[i] == ',':
            if not have_value or comma_pending:
                raise records.error(section, line, "explicit values; null values are unsupported")
            comma_pending = True
            i += 1
            continue
        start = i
        repeat = 1
        match = re.match(r"([0-9]+)\*", line[i:])
        if match:
            try:
                repeat = int(match[1])
            except ValueError as exc:
                raise records.error(section, match[0], "a representable positive repetition count") from exc
            i += len(match[0])
            if repeat == 0 or i == len(line) or line[i].isspace() or line[i] in ',/':
                raise records.error(section, line[start:], "positive repetition count and explicit value (null repetition unsupported)")
        quoted = line[i] in "\"'"
        if quoted:
            quote = line[i]
            i += 1
            chars = []
            while i < len(line):
                if line[i] == quote:
                    if i + 1 < len(line) and line[i + 1] == quote:
                        chars.append(quote)
                        i += 2
                        continue
                    i += 1
                    break
                chars.append(line[i])
                i += 1
            else:
                raise records.error(section, line[start:], "a closed quoted character value on one record")
            value = ''.join(chars)
            if i < len(line) and not line[i].isspace() and line[i] != ',':
                raise records.error(section, line[start:], "a separator after a quoted value")
        else:
            while i < len(line) and not line[i].isspace() and line[i] != ',':
                i += 1
            value = line[start + (len(match[0]) if match else 0):i]
            if any(c in value for c in "/!*()\"'"):
                raise records.error(section, line[start:i], "a scalar value; slash, comments, nulls and complex syntax are unsupported")
        have_value = True
        comma_pending = False
        yield _Token(line[start:i], value, records.position, quoted, repeat)


class _Read:
    def __init__(self, records: _Records, section: str) -> None:
        self.records = records
        self.section = section
        self.tokens: Iterator[_Token] = iter(())
        self.pending: _Token | None = None
        self.remaining = 0
        self.last: _Token | None = None
        self.fields: dict[str, _Token] = {}

    def token(self) -> _Token:
        if self.remaining:
            self.remaining -= 1
            return self.pending  # type: ignore[return-value]
        while True:
            token = next(self.tokens, None)
            if token is not None:
                self.pending = token
                self.remaining = token.repeat - 1
                self.last = token
                return token
            line = self.records.next(self.section)
            self.tokens = iter(tuple(_lex(self.records, self.section, line)))

    def fail(self, token: _Token, expected: str) -> LegacyInputParseError:
        return self.records.error(self.section, token.raw, expected, token.record)

    def integer(self, name: str, *, minimum: int | None = None,
                maximum: int | None = None, allowed: tuple[int, ...] | None = None) -> int:
        token = self.token()
        self.fields[name] = token
        value_text = token.value
        expected = f"{name}: an integer"
        if minimum is not None and maximum is not None:
            expected += f" in {minimum}..{maximum}"
        elif minimum is not None:
            expected += f" >= {minimum}"
        if allowed is not None:
            expected += f" in {allowed}"
        if token.quoted or not _INT_RE.fullmatch(value_text):
            raise self.fail(token, expected)
        try:
            value = int(value_text)
        except ValueError as exc:
            raise self.fail(token, expected) from exc
        if ((minimum is not None and value < minimum)
                or (maximum is not None and value > maximum)
                or (allowed is not None and value not in allowed)):
            raise self.fail(token, expected)
        return value

    def real(self, name: str, *, positive: bool = False) -> float:
        token = self.token()
        self.fields[name] = token
        value_text = token.value
        expected = f"{name}: a finite {'positive ' if positive else ''}real number"
        if token.quoted or not re.fullmatch(
            r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[EeDd][+-]?[0-9]+|[+-][0-9]+)?",
            value_text,
        ):
            raise self.fail(token, expected)
        value_text = re.sub(r"[Dd]", "e", value_text)
        # Fortran's exponent-letter-omitted form, e.g. 1.25-2.
        if not re.search('[Ee]', value_text):
            value_text = re.sub(r"(?<=[0-9.])([+-][0-9]+)$", r"e\1", value_text)
        value = float(value_text)
        if not math.isfinite(value) or (positive and value <= 0):
            raise self.fail(token, expected)
        return value

    def flag(self, name: str, allowed: tuple[int, ...] = (0, 1)) -> int:
        return self.integer(name, allowed=allowed)

    def model(self, factory: Callable[..., T], **kwargs: object) -> T:
        try:
            return factory(**kwargs)
        except (ValueError, TypeError, OverflowError) as exc:
            token = next((token for name, token in self.fields.items()
                          if name in str(exc)), self.last)
            assert token is not None
            line = self.records.lines[token.record - 1]
            raise self.records.error(self.section, line, str(exc), token.record) from exc


def parse_for007(path: str | Path) -> LegacyInput:
    """Parse through EXPINT using one cursor per current source READ.

    Later scan/refinement sections remain outside this first-stage parser.
    Formatted character records are never passed through the list lexer.
    """
    records = _Records(Path(path))
    title_raw = records.next("title")
    read = records.read("beam voltage and zone axis")
    energy = read.real("beam voltage", positive=True)
    refined = read.flag("energy refinement flag")
    zone = tuple(read.integer("zone-axis index") for _ in range(3))
    read = records.read("slice count and mesh dimensions")
    count = read.integer("slice count", minimum=1)
    mesh = read.model(MeshConfig, mesh_x=read.integer("mesh_x", minimum=1),
                      mesh_y=read.integer("mesh_y", minimum=1))
    slices = []
    for i in range(count):
        read = records.read(f"slice geometry {i + 1}")
        slices.append(read.model(
            SliceGeometry,
            cell_a_angstrom=read.real("cell_a_angstrom", positive=True),
            cell_b_angstrom=read.real("cell_b_angstrom", positive=True),
            lattice_gamma_degrees=read.real("lattice_gamma_degrees"),
            thickness_angstrom=read.real("thickness_angstrom", positive=True),
            x_shift_fractional=read.real("x_shift_fractional"),
            y_shift_fractional=read.real("y_shift_fractional"),
        ))
    species_count = records.read("species count").integer("species count", minimum=1)
    species = []
    for i in range(species_count):
        raw = records.next(f"species {i + 1} label")
        # FORMAT(A2), not list-directed: quotes and spaces are literal.
        label = raw[:2].ljust(2)
        read = records.read(f"species {i + 1} atomic number")
        species.append(read.model(Species, label=label, label_raw=raw,
                                  atomic_number=read.integer("atomic number", minimum=1)))
    records.read("DSAP").real("DSAP")
    records.read("DSMAX").real("DSMAX")
    regions = records.read("region count").integer("region count", minimum=1)
    read = records.read("region depths and masks")
    for _ in range(regions):
        read.real("region depth", positive=True)
        read.flag("thickness mask")
    read = records.read("region slice counts")
    region_slices = [read.integer("region slice count", minimum=1) for _ in range(regions)]
    if sum(region_slices) != count:
        raise read.fail(read.last, f"region slice counts summing to slice count {count}")
    for i in range(regions):
        flag = records.read(f"region {i + 1} vacancy flag").flag("vacancy flag", (0, 1, 2))
        read = records.read(f"region {i + 1} vacancy parameters")
        for _ in range(2 if flag else 1):
            read.real("vacancy coefficient")
            read.flag("vacancy refinement mask")
    read = records.read("lattice dilation")
    read.real("dilation")
    read.flag("dilation refinement mask")
    absorption_raw = records.next("absorption flag")
    if absorption_raw[:2] == "YE":
        for i in range(species_count):
            read = records.read(f"absorption coefficients {i + 1}")
            for coefficient in ('A', 'B', 'C'):
                read.real(f"absorption {coefficient}")
                read.flag(f"absorption {coefficient} mask", (-1, 0, 1) if coefficient == 'C' else (0, 1))
    for i in range(regions):
        read = records.read(f"Debye-Waller factors region {i + 1}")
        for _ in range(species_count):
            read.real("Debye-Waller factor")
        read = records.read(f"Debye-Waller masks region {i + 1}")
        for _ in range(species_count):
            read.flag("Debye-Waller mask")
        read = records.read(f"free-surface parameters region {i + 1}")
        read.flag("free-surface type", (0, 1, 2))
        for _ in range(3):
            read.real("free-surface coefficient")
    atoms_by_slice = []
    for i in range(count):
        read = records.read(f"slice {i + 1} atom header")
        atoms_count = read.integer("atom count", minimum=0)
        kind = read.token()
        if kind.value not in ('S', 'G'):
            raise read.fail(kind, "record type exactly S or G")
        atoms = []
        for j in range(atoms_count):
            read = records.read(f"slice {i + 1} atom {j + 1}")
            index = read.integer("species index", minimum=1, maximum=species_count)
            x = read.real("x_fractional")
            x_mask = read.flag("x mask") if kind.value == 'S' else 0
            y = read.real("y_fractional")
            y_mask = read.flag("y mask") if kind.value == 'S' else 0
            occupancy = read.real("starting_occupancy")
            atoms.append(read.model(AtomSite, species_index_one_based=index,
                                    x_fractional=x, y_fractional=y,
                                    starting_occupancy=occupancy, record_type=kind.value,
                                    record_type_raw=kind.raw, x_mask=x_mask, y_mask=y_mask))
        atoms_by_slice.append(tuple(atoms))
    output_raw = records.next("output filename BMNAME")
    if not output_raw.strip():
        raise records.error("output filename BMNAME", output_raw, "a non-empty filename")
    # A single READ consumes MBOUT and its dependent implied-DO list.
    read = records.read("output beam count and indices")
    beams = read.integer("output beam count", minimum=0)
    for _ in range(2 * beams):
        read.integer("beam index")
    experimental_raw = records.next("experimental intensity filename EXPINT")
    if not experimental_raw.strip():
        raise records.error("experimental intensity filename EXPINT", experimental_raw, "a non-empty filename")
    return LegacyInput(
        title=title_raw.rstrip(), energy_kv=energy, energy_refined=refined == 1,
        zone_axis=zone, mesh=mesh, slices=tuple(slices), species=tuple(species),
        atoms_by_slice=tuple(atoms_by_slice), output_filename=output_raw.rstrip(),
        experimental_intensity_filename=experimental_raw.rstrip(), title_raw=title_raw,
        absorption_raw=absorption_raw, output_filename_raw=output_raw,
        experimental_intensity_filename_raw=experimental_raw,
    )

"""Initial legacy scan records and explicit undilated reciprocal-cell tilts."""

from dataclasses import dataclass
import math
from numbers import Integral, Real
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .legacy_io import _Records, _parse_for007_records


def _coordinates(values, name, size, nonzero_xy=False):
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must contain {size} finite real values")
    try:
        result = tuple(values)
    except TypeError as exc:
        raise TypeError(f"{name} must contain {size} finite real values") from exc
    if len(result) != size:
        raise ValueError(f"{name} must contain {size} finite real values")
    for value in result:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
            raise TypeError(f"{name} must contain finite real values excluding bool")
        if not math.isfinite(value):
            raise ValueError(f"{name} must contain finite real values")
    result = tuple(float(value) for value in result)
    if nonzero_xy and result[:2] == (0., 0.):
        raise ValueError(f"{name} must have a nonzero in-plane component; source normalization would be singular")
    return result


@dataclass(frozen=True)
class LegacyScanDefinition:
    """Initial scan-basis corners and counts, with no refinement or dilation.

    basis_1=(AY,BEE,NOTHING), basis_2=(CEE,DEE,NAUGHT). Third components
    are retained but unused by the readable source transformation. Corners
    contain dimensionless coefficients of those two scan basis vectors.
    scan_columns=NUMP, scan_rows=NUML; neither is a wave mesh dimension.
    """

    basis_1: tuple[float, float, float]
    basis_2: tuple[float, float, float]
    corner_a: tuple[float, float]
    corner_b: tuple[float, float]
    corner_c: tuple[float, float]
    scan_columns: int
    scan_rows: int

    def __post_init__(self) -> None:
        for name in ("basis_1", "basis_2", "corner_a", "corner_b", "corner_c"):
            basis = name.startswith("basis")
            object.__setattr__(self, name, _coordinates(getattr(self, name), name,
                                                       3 if basis else 2, basis))
        for name in ("scan_columns", "scan_rows"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
                raise TypeError(f"{name} must be a positive integer excluding bool")
            if value <= 0:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, int(value))


class _HeaderRecords(_Records):
    """Stream just one list-directed header READ, never its intensity body."""

    def __init__(self, filename, stream):
        self.filename, self.stream = filename, stream
        self.position, self.lines = 0, []

    def next(self, section):
        line = self.stream.readline()
        if not line:
            raise self.error(section, "<EOF>", "another input record", self.position+1)
        try:
            line = line.decode("utf-8")
        except UnicodeError as exc:
            raise self.error(section, repr(line), "a UTF-8 header record", self.position+1) from exc
        self.position += 1
        self.lines.append(line.rstrip("\r\n"))
        return self.lines[-1]


def parse_for007_scan(
    path: str | Path, *, experimental_header_path: str | Path | None = None,
) -> LegacyScanDefinition:
    """Read initial basis/corners plus NUMP,NUML from the named data header.

    The default header path is EXPINT relative to the FOR007 directory, an
    explicit Python file-location policy. An override is relative to the
    caller's working directory. Only NUMP,NUML,GRADIENTOPT are consumed there;
    the third integer is discarded without interpreting its objective option.
    FOR007 parsing stops after C, before MASKT or any later refinement records.
    """
    records = _Records(Path(path))
    prefix = _parse_for007_records(records)
    header_path = (Path(experimental_header_path) if experimental_header_path is not None
                   else Path(path).parent / prefix.experimental_intensity_filename)
    # Binary lines avoid decoding any buffered experimental-intensity body.
    try:
        with header_path.open("rb") as stream:
            header = _HeaderRecords(header_path, stream)
            read = header.read("scan dimensions header")
            columns = read.integer("scan_columns", minimum=1)
            rows = read.integer("scan_rows", minimum=1)
            read.integer("GRADIENTOPT")
    except (OSError, UnicodeError) as exc:
        raise records.error("scan dimensions header", str(header_path),
                            "a readable header containing NUMP,NUML,GRADIENTOPT") from exc
    read = records.read("scan basis")
    basis = []
    for name in ("basis_1", "basis_2"):
        values = tuple(read.real(f"{name}[{i}]") for i in range(3))
        basis.append(read.model(_coordinates, values=values, name=name, size=3, nonzero_xy=True))
    corners = []
    for name in ("corner_a", "corner_b", "corner_c"):
        read = records.read(f"scan {name}")
        corners.append(tuple(read.real(f"{name}[{i}]") for i in range(2)))
    return read.model(LegacyScanDefinition, basis_1=basis[0], basis_2=basis[1],
        corner_a=corners[0], corner_b=corners[1], corner_c=corners[2],
        scan_columns=columns, scan_rows=rows)


def generate_scan_tilts(scan: LegacyScanDefinition) -> NDArray[np.float64]:
    """Return (scan_rows*scan_columns,2) float64, columns H,K; columns fastest.

    Algebraically cancel NORMFACTOR/INVNORMFACTOR and leave dilation to the
    existing forward API. No singleton rule is invented for source 0/0.
    """
    if not isinstance(scan, LegacyScanDefinition):
        raise TypeError("scan must be LegacyScanDefinition")
    if scan.scan_columns == 1 or scan.scan_rows == 1:
        raise ValueError("singleton scan dimensions are undefined in the source: NUMP-1 and NUML-1 denominators; both must exceed one")
    a, b, c = (np.asarray(value, dtype=np.float64) for value in
               (scan.corner_a, scan.corner_b, scan.corner_c))
    row = np.arange(scan.scan_rows, dtype=np.float64)[:, None]
    col = np.arange(scan.scan_columns, dtype=np.float64)[None, :]
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            step_row, step_col = (c-a)/(scan.scan_rows-1), (b-a)/(scan.scan_columns-1)
            u = a[0] + row*step_row[0] + col*step_col[0]
            v = a[1] + col*step_col[1] + row*step_row[1]
            h = scan.basis_1[0]*u + scan.basis_2[0]*v
            k = scan.basis_1[1]*u + scan.basis_2[1]*v
            result = np.stack((h, k), axis=-1).reshape(-1, 2)
    except (FloatingPointError, OverflowError) as exc:
        raise ValueError("generated scan tilts must remain finite in float64") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError("generated scan tilts must remain finite in float64")
    return result

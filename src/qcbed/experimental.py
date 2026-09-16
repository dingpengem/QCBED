"""Unmodified experimental file values in public (tilt, beam) order."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .legacy_io import _Records
from .multislice import BeamIndices, _beams, _integer, _snapshot
from .scan import LegacyScanDefinition


@dataclass(frozen=True, eq=False)
class ExperimentalIntensityData:
    """Raw XPINTS, VARIANCE and WEIGHTS; no scaling or weighting applied.

    All arrays have float64 shape (scan_rows*scan_columns, beam_count).
    beam_indices_hk must be supplied in the external FOR007 MBOUT list order.
    gradient_option preserves the header integer without running GRADIENT.
    """

    scan_columns: int
    scan_rows: int
    gradient_option: int
    beam_indices_hk: BeamIndices
    intensities: NDArray[np.float64]
    variances: NDArray[np.float64]
    weights: NDArray[np.float64]

    def __post_init__(self) -> None:
        for name in ("scan_columns", "scan_rows", "gradient_option"):
            value = _integer(getattr(self, name), name)
            if name != "gradient_option" and value <= 0:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, value)
        beams = _beams(self.beam_indices_hk)
        object.__setattr__(self, "beam_indices_hk", beams)
        shape = (self.scan_rows*self.scan_columns, len(beams))
        for name in ("intensities", "variances", "weights"):
            values = getattr(self, name)
            if not isinstance(values, np.ndarray) or values.dtype != np.float64:
                raise TypeError(f"{name} must be a float64 NumPy array")
            if values.shape != shape:
                raise ValueError(f"{name} must have (tilt, beam) shape {shape}")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{name} must contain finite values")
            object.__setattr__(self, name, _snapshot(values))


def parse_experimental_intensities(
    path: str | Path, *, beam_indices_hk: BeamIndices,
    scan: LegacyScanDefinition | None = None,
) -> ExperimentalIntensityData:
    """Parse the unit-12 READ sequence, preserving each logical row READ.

    Beam count/order are absent from this file and must be supplied explicitly.
    Optional scan validates both dimensions without generating or deduplicating
    tilts. Each block consists of beam_count*NUML READs of NUMP real values.
    Advancing READs permit continuation and discard surplus values on their
    final record. Extra nonblank records after all three blocks are rejected.
    """
    beams = _beams(beam_indices_hk)
    if scan is not None and not isinstance(scan, LegacyScanDefinition):
        raise TypeError("scan must be LegacyScanDefinition or None")
    records = _Records(Path(path))
    read = records.read("experimental header")
    columns = read.integer("scan_columns NUMP", minimum=1)
    rows = read.integer("scan_rows NUML", minimum=1)
    option = read.integer("gradient_option GRADIENTOPT")
    if scan is not None and (columns, rows) != (scan.scan_columns, scan.scan_rows):
        raise read.fail(read.last, f"scan dimensions NUMP,NUML = {scan.scan_columns},{scan.scan_rows}; "
                        f"observed {columns},{rows}")
    blocks = {}
    for name in ("intensities", "variances", "weights"):
        block = np.empty((rows*columns, len(beams)), dtype=np.float64)
        for beam in range(len(beams)):
            for row in range(rows):
                read = records.read(f"experimental {name}, beam {beam+1}, scan row {row+1}")
                for col in range(columns):
                    block[row*columns+col, beam] = read.real(f"{name} column {col+1}")
        blocks[name] = block
    for index in range(records.position, len(records.lines)):
        if records.lines[index].strip():
            raise records.error("experimental trailing records", records.lines[index],
                                "end of data after intensities, variances and weights", index+1)
    return ExperimentalIntensityData(columns, rows, option, beams, **blocks)

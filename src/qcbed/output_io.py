"""Source-derived unit-6 output: scaled objective-space values, no header."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .legacy_io import _Records, _lex
from .multislice import _beams, _snapshot, _integer
from .scan import LegacyScanDefinition, generate_scan_tilts
from .intensity_scaling import _scalar, apply_legacy_intensity_scaling


def _array(values, shape, name):
    if not isinstance(values,np.ndarray) or values.dtype != np.float64:
        raise TypeError(f"{name} must be a float64 array")
    if values.shape != shape or not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return _snapshot(values)


def reshape_scan_map(values, scan: LegacyScanDefinition):
    """Row r, column c receives tilt r*NUMP+c; no transpose or display scaling."""
    values = _array(values,(scan.scan_rows*scan.scan_columns,),"values")
    return _snapshot(values.reshape(scan.scan_rows,scan.scan_columns))


@dataclass(frozen=True, eq=False)
class LegacySimulatedOutput:
    """Headerless unit-6 values; scan and beam identities supplied externally.

    Raw intensity, GRADIENT option and YNORM cannot be recovered from this file.
    """
    scan: LegacyScanDefinition
    beam_indices_hk: tuple
    scaled_intensity: np.ndarray

    def __post_init__(self):
        if not isinstance(self.scan,LegacyScanDefinition):
            raise TypeError("scan must be LegacyScanDefinition")
        beams = _beams(self.beam_indices_hk)
        object.__setattr__(self,"beam_indices_hk",beams)
        object.__setattr__(self,"scaled_intensity",_array(self.scaled_intensity,
            (self.scan.scan_rows*self.scan.scan_columns,len(beams)),"scaled_intensity"))


@dataclass(frozen=True, eq=False)
class SimulatedCbedData:
    """Raw |amplitude|², post-GRADIENT values and globally YNORM-scaled output.

    Arrays use (tilt, beam); tilts are undilated fractional reciprocal H,K.
    """
    scan: LegacyScanDefinition
    beam_indices_hk: tuple
    raw_intensity: np.ndarray
    objective_intensity: np.ndarray
    ynorm: float
    gradient_option: int
    objective: float
    scaled_intensity: np.ndarray = field(init=False)
    tilts: np.ndarray = field(init=False)

    def __post_init__(self):
        if not isinstance(self.scan,LegacyScanDefinition):
            raise TypeError("scan must be LegacyScanDefinition")
        beams = _beams(self.beam_indices_hk)
        object.__setattr__(self,"beam_indices_hk",beams)
        shape=(self.scan.scan_rows*self.scan.scan_columns,len(beams))
        for name in ("raw_intensity","objective_intensity"):
            object.__setattr__(self,name,_array(getattr(self,name),shape,name))
        option = _integer(self.gradient_option,"gradient_option")
        if option not in (0,1):
            raise ValueError("gradient_option must be 0 or 1")
        # The source objective may be undefined for zero experimental variance.
        # Keep that diagnostic; finite maps/YNORM remain mandatory for output.
        for name in ("ynorm","objective"):
            object.__setattr__(self,name,_scalar(getattr(self,name),name,finite=name=="ynorm"))
        object.__setattr__(self,"scaled_intensity",_array(
            apply_legacy_intensity_scaling(self.objective_intensity,self.ynorm),shape,"scaled_intensity"))
        object.__setattr__(self,"tilts",_snapshot(generate_scan_tilts(self.scan)))


def write_legacy_out_dat(path, data: SimulatedCbedData | LegacySimulatedOutput):
    """Beam/row/column order with blank row separators; portable .17E precision.

    This is numeric source-format compatibility, not compiler byte equivalence.
    """
    with Path(path).open("w",encoding="ascii",newline="\n") as stream:
        for b in range(len(data.beam_indices_hk)):
            for row in reshape_scan_map(data.scaled_intensity[:,b],data.scan):
                stream.write(" ".join(f"{value:.17E}" for value in row)+"\n\n")


def parse_legacy_out_dat(path, *, scan, beam_indices_hk):
    """Strict single dataset, accepting runtime line wrapping and blank lines.

    Numeric tokens are counted across records because WRITE(*), not READ row
    boundaries, defines this headerless format. No surplus snapshot is discarded.
    """
    beams=_beams(beam_indices_hk)
    records=_Records(Path(path))
    count=scan.scan_rows*scan.scan_columns*len(beams)
    # Validate/count the numeric stream before reading: refuse extra snapshots.
    total=0
    for number,line in enumerate(records.lines,1):
        records.position=number
        total+=sum(token.repeat for token in _lex(records,"simulated output",line))
    if total!=count:
        raise records.error("simulated output",f"{total} values",f"exactly {count} values (one dataset)")
    records.position=0
    read=records.read("simulated output")
    flat=np.array([read.real("scaled objective-space value") for _ in range(count)],dtype=np.float64)
    values=flat.reshape(len(beams),-1).T
    return LegacySimulatedOutput(scan,beams,values)

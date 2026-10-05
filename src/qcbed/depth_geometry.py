"""Canonical strict-< whole-slice depth expansion, preserving addition order."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class RequestedDepthGeometry:
    requested_depth: float
    effective_depth: float
    sequence: tuple[int, ...]
    stored_slice_count: int
    plateau_lower_exclusive: float

    @property
    def occurrence_count(self):return len(self.sequence)


def requested_depth_to_sample_geometry(requested_depth, slice_thicknesses):
    """Cycle stored slices while accumulated depth < request; never interpolate.

    Accumulate sequentially, not ceil(request/thickness) or count*thickness:
    their floating-point boundaries need not equal the source loop's boundaries.
    The current sequence owns (previous accumulated depth, effective_depth].
    """
    if isinstance(requested_depth,bool) or not math.isfinite(requested_depth) or requested_depth<=0:
        raise ValueError('requested depth must be finite and positive')
    thicknesses=tuple(slice_thicknesses)
    if not thicknesses or any(isinstance(t,bool) or not math.isfinite(t) or t<=0 for t in thicknesses):
        raise ValueError('slice thicknesses must be finite and positive')
    sequence=[]; total=0.; previous=0.
    while total<requested_depth:
        i=len(sequence)%len(thicknesses)
        previous=total
        updated=total+thicknesses[i]
        if updated<=total:raise ValueError('slice thickness must advance accumulated depth')
        sequence.append(i);total=updated
    return RequestedDepthGeometry(float(requested_depth),total,tuple(sequence),len(thicknesses),previous)

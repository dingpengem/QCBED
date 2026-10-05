"""Deterministic direct-cell replication shared by 3-D viewing and projection."""
from itertools import product
from numbers import Integral
import numpy as np
from .lattice import _readonly


def replica_offsets(repeats):
    """Integer cell offsets in lexicographic order, c fastest; no centering/wrapping."""
    try:
        repeats = tuple(repeats)
    except TypeError as exc:
        raise ValueError('repeats must contain three positive integers') from exc
    if len(repeats) != 3 or any(isinstance(v, (bool, np.bool_)) or
            not isinstance(v, Integral) or v <= 0 for v in repeats):
        raise ValueError('repeats must contain three positive integers')
    return tuple(product(*(range(int(n)) for n in repeats)))


def replicated_cartesian_positions(structure, offsets):
    """Cell outermost, original site innermost; output Å, with no site merging."""
    fractional = np.array([s.fractional_position for s in structure.sites])
    return _readonly((fractional[None,:,:]+np.asarray(offsets)[:,None,:]).reshape(-1,3)
                     @ structure.direct_lattice.vectors)

"""Explicit [001] legacy compatibility coordinates, never canonical indexing.

The legacy reciprocal directions are conventional (110) and (1,-1,0).
All transformations are exact integer algebra. No beam filtering is performed.
"""
from .orientation import _indices


def legacy_for007_001_basis_to_conventional_hkl(p, q):
    """Convert legacy integer (p,q) to conventional (p+q,p-q,0)."""
    p, q, _ = _indices((p, q, 0), 'legacy (p,q)')
    return p + q, p - q, 0


def conventional_hkl_to_legacy_for007_001_basis(hkl):
    """Reject non-ZOLZ or parity-incompatible reflections; never round."""
    h, k, l = _indices(hkl, 'conventional hkl')
    if l != 0 or (h + k) % 2:
        raise ValueError('legacy [001] requires l=0 and same-parity h,k')
    return (h + k) // 2, (h - k) // 2


def fcc_conventional_reflection_diagnostic(hkl):
    """Return direct/allowed/forbidden, for reporting only, not dynamics."""
    h, k, l = _indices(hkl, 'conventional hkl')
    if (h, k, l) == (0, 0, 0):
        return 'direct'
    return 'allowed' if h % 2 == k % 2 == l % 2 else 'forbidden'

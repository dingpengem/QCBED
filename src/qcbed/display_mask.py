"""Experimental display masks and counts; never used by objective arithmetic."""

import numpy as np

from .multislice import _snapshot


def valid_display_mask(experimental):
    """Immutable (tilt, beam) boolean mask: every nonzero weight is visible."""
    return _snapshot(experimental.weights != 0)


def experimental_mask_statistics(experimental):
    """Counts from file values, including any visible zero-variance locations.

    Indices in location records are zero-based; beam_hk retains legacy H/K.
    Fractional or negative nonzero weights are not rounded into binary weights.
    """
    e = experimental
    zero = e.variances == 0
    masked = e.weights == 0

    def counts(z, m):
        return dict(n_total=int(z.size), n_zero_variance=int(z.sum()),
                    n_zero_weight=int(m.sum()),
                    n_zero_variance_and_zero_weight=int((z & m).sum()),
                    n_zero_variance_and_nonzero_weight=int((z & ~m).sum()),
                    n_nonzero_variance_and_zero_weight=int((~z & m).sum()))

    result = counts(zero, masked)
    result.update(unique_weights=np.unique(e.weights).tolist(),
                  binary_weights=bool(np.all((e.weights == 0) | (e.weights == 1))),
                  n_negative_weights=int(np.count_nonzero(e.weights < 0)),
                  n_nonfinite_weights=int(np.count_nonzero(~np.isfinite(e.weights))),
                  same_mask_all_beams=bool(np.all(masked == masked[:, :1])),
                  per_beam=[dict(beam_hk=beam, **counts(zero[:, b], masked[:, b]))
                            for b, beam in enumerate(e.beam_indices_hk)],
                  visible_zero_variance=[dict(beam_index=int(b), beam_hk=e.beam_indices_hk[b],
                      tilt_index=int(t), row=int(t // e.scan_columns), column=int(t % e.scan_columns),
                      intensity=float(e.intensities[t,b]), variance=float(e.variances[t,b]),
                      weight=float(e.weights[t,b])) for t,b in np.argwhere(zero & ~masked)])
    return result

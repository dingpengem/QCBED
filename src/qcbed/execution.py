"""Immutable forward-execution policy, independent of objective and physics."""

from dataclasses import dataclass
from typing import Literal

from .tilt_series import _batch_size


@dataclass(frozen=True)
class ForwardExecutionConfig:
    """Serial reference by default; batched NumPy uses 64 tilts by default.

    Process execution remains available through the standalone tilt API only.
    Serial rejects batch settings rather than silently ignoring them.
    """

    execution_backend: Literal["serial", "batched"] = "serial"
    batch_size: int | None = None

    def __post_init__(self):
        if not isinstance(self.execution_backend, str) or self.execution_backend not in ("serial", "batched"):
            raise ValueError("execution_backend must be 'serial' or 'batched'")
        if self.execution_backend == "serial":
            if self.batch_size is not None:
                raise ValueError("serial execution requires batch_size=None")
        else:
            object.__setattr__(self, "batch_size", _batch_size(
                64 if self.batch_size is None else self.batch_size))


def _execution(value):
    if value is None:
        return ForwardExecutionConfig()
    if not isinstance(value, ForwardExecutionConfig):
        raise TypeError("execution must be ForwardExecutionConfig")
    return value

"""Optional maintainer-only evidence kept entirely outside the shared checkout."""
import os
from pathlib import Path
import pytest


def private_reference_root():
    value=os.environ.get('QCBED_PRIVATE_REFERENCE_ROOT')
    if not value:
        pytest.skip('Optional private-reference evidence unavailable; set QCBED_PRIVATE_REFERENCE_ROOT')
    root=Path(value).resolve()
    assert root.is_dir(), 'Configured private-reference root is not a directory'
    assert not root.is_relative_to(Path(__file__).resolve().parents[1]), 'Private evidence must remain outside the shared checkout'
    return root

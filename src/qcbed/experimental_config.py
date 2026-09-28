"""Explicit configuration-driven experimental I/O, separate from compilation."""
from pathlib import Path

from .config_compiler import CompiledConfig
from .experimental import ExperimentalIntensityData, _parse_experimental_body
from .legacy_io import _Records


def load_experimental_data(
    config: CompiledConfig, *, base_directory: str | Path,
) -> ExperimentalIntensityData:
    """Read configured body blocks, optionally validating a legacy logical header.

    Relative paths are relative to the explicitly supplied base directory.
    No header autodetection, metadata overrides, preprocessing or normalization.
    The same advancing READ implementation is used by the legacy reader.
    """
    if not isinstance(config, CompiledConfig):
        raise TypeError('config must be CompiledConfig')
    spec = config.experiment
    columns, rows = config.state.scan.scan_columns, config.state.scan.scan_rows
    option = spec.gradient_option
    path = Path(spec.path)
    if not path.is_absolute():
        path = Path(base_directory)/path
    records = _Records(path)
    if spec.header == 'validate':
        read = records.read('experimental header')
        observed = (read.integer('scan_columns NUMP', minimum=1),
                    read.integer('scan_rows NUML', minimum=1),
                    read.integer('gradient_option GRADIENTOPT'))
        for name, actual, expected in zip(
                ('scan_columns NUMP', 'scan_rows NUML', 'gradient_option GRADIENTOPT'),
                observed, (columns, rows, option)):
            if actual != expected:
                raise read.fail(read.fields[name], f'{name} = {expected} from canonical configuration; observed {actual}')
    return _parse_experimental_body(records, columns, rows, option, spec.beam_indices_hk)

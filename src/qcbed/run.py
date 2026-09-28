"""Canonical production entry point; delegates all refinement to full_runner."""
import argparse
import json
from pathlib import Path

from .yaml_config import load_yaml_config


def run_yaml_refinement(path: str | Path, *, resume=False, progress_hook=None):
    """Run exactly the YAML settings; resume only matching canonical checkpoints.

    Identity uses the resolved compiled configuration and parsed experimental
    arrays, never YAML bytes or FOR007. No quick-mode or implicit overrides.
    """
    from .canonical_run import run_config
    return run_config(load_yaml_config(path), base_directory=Path(path).resolve().parent,
                      resume=resume, progress_hook=progress_hook)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Run canonical QCBED refinement without FOR007 input.')
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--resume',action='store_true',help='Resume the configured output directory.')
    args = parser.parse_args(argv)
    result = run_yaml_refinement(args.config,resume=args.resume)
    print(json.dumps(result.summary,indent=2))
    return 0 if result.status == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())

"""Explicit compatibility CLI; canonical production runs never import this module."""
import argparse
import json
from pathlib import Path

from .yaml_config import load_yaml_config
from .legacy_config import MigrationPolicy, import_legacy_config, write_yaml_config, export_for007, _plain
from .configuration import RunConfig, BoundsConfig
from .execution import ForwardExecutionConfig


def main(argv=None):
    parser = argparse.ArgumentParser(description='One-time legacy migration and optional FOR007 export')
    sub = parser.add_subparsers(dest='command',required=True)
    export = sub.add_parser('export')
    export.add_argument('--config',required=True)
    export.add_argument('--destination',required=True)
    export.add_argument('--bundle',action='store_true')
    export.add_argument('--experimental-filename')
    export.add_argument('--output-filename')
    export.add_argument('--filename-limits',nargs=2,type=int,metavar=('EXPINT','BMNAME'))
    imp = sub.add_parser('import')
    imp.add_argument('--for007',required=True)
    imp.add_argument('--experimental',required=True)
    imp.add_argument('--yaml',required=True)
    imp.add_argument('--seed',type=int,required=True)
    imp.add_argument('--execution',choices=['serial','batched'],required=True)
    imp.add_argument('--batch-size',type=int)
    imp.add_argument('--objective',choices=['legacy','masked'],required=True)
    imp.add_argument('--invalid-trial-policy',choices=['raise','reject'],required=True)
    imp.add_argument('--bounds',choices=['none'],required=True,help='Explicitly no imported hard bounds; add bounds in YAML after migration.')
    imp.add_argument('--output-directory',required=True)
    imp.add_argument('--output-filename',required=True)
    imp.add_argument('--snapshot-every',type=int,required=True,help='0 disables snapshots')
    imp.add_argument('--checkpoint-every',type=int,required=True)
    imp.add_argument('--max-evaluations',type=int,required=True,help='0 means unlimited')
    imp.add_argument('--max-seconds',type=float,required=True,help='0 means unlimited')
    args = parser.parse_args(argv)
    if args.command == 'export':
        report = export_for007(load_yaml_config(args.config),args.destination,bundle=args.bundle,
            experimental_filename=args.experimental_filename,output_filename=args.output_filename,filename_limits=args.filename_limits)
        print(json.dumps(_plain(report),indent=2))
    else:
        policy = MigrationPolicy(args.seed,ForwardExecutionConfig(args.execution,args.batch_size),
            RunConfig(str(Path(args.output_directory).resolve()),str(Path(args.output_filename).resolve()),
                args.snapshot_every or None,args.checkpoint_every,args.max_evaluations or None,args.max_seconds or None),
            args.objective,args.invalid_trial_policy,BoundsConfig())
        raw = import_legacy_config(args.for007,experimental_path=args.experimental,policy=policy)
        write_yaml_config(raw,args.yaml)
        print(f'Wrote {args.yaml}; ordered IDs generated; Python policies explicitly supplied. Experimental mode: validate.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

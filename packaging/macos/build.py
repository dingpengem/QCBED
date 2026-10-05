"""Reproducible-location clean build; leaves the development environment intact."""
import json
import hashlib
from datetime import datetime, timezone
import os
from pathlib import Path
import platform
import subprocess
import sys
from importlib.metadata import version


def main():
    if sys.platform != 'darwin' or platform.machine() != 'arm64':
        raise SystemExit('Use native Apple Silicon macOS Python')
    root = Path(__file__).resolve().parents[2]
    build = root/'build/macos'
    build.mkdir(parents=True,exist_ok=True)
    env = dict(os.environ, PYINSTALLER_CONFIG_DIR=str(build/'cache'), MPLCONFIGDIR=str(build/'mpl-cache'))
    subprocess.run([sys.executable,'-m','PyInstaller','--clean','--noconfirm',
        '--workpath',str(build/'work'),'--distpath',str(root/'dist'),
        str(root/'packaging/macos/QCBED.spec')],env=env,check=True)
    source_hash=hashlib.sha256()
    for path in sorted([*(root/'src/qcbed').rglob('*.py'),*(root/'packaging/macos').glob('*.py'),root/'packaging/macos/QCBED.spec']):
        source_hash.update(str(path.relative_to(root)).encode());source_hash.update(path.read_bytes())
    names = ('qcbed','pyinstaller','pyinstaller-hooks-contrib','PySide6','numpy','scipy',
             'h5py','matplotlib','psutil','PyYAML','ase')
    (build/'build-environment.json').write_text(json.dumps(dict(
        source_sha256=source_hash.hexdigest(),built_at_utc=datetime.now(timezone.utc).isoformat(),
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        git_dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True).strip()),
        python=sys.version,executable=sys.executable,platform=platform.platform(),
        architecture=platform.machine(),versions={name:version(name) for name in names}),indent=2)+'\n')
    print(root/'dist/QCBED.app')


if __name__ == '__main__': main()

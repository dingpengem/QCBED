"""Reproducible-location clean build; leaves the development environment intact."""
import json
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
    names = ('qcbed','pyinstaller','pyinstaller-hooks-contrib','PySide6','numpy','scipy',
             'h5py','matplotlib','psutil','PyYAML')
    (build/'build-environment.json').write_text(json.dumps(dict(
        python=sys.version,executable=sys.executable,platform=platform.platform(),
        architecture=platform.machine(),versions={name:version(name) for name in names}),indent=2)+'\n')
    print(root/'dist/QCBED.app')


if __name__ == '__main__': main()

# Build on macOS arm64: python packaging/macos/build.py
from pathlib import Path
import platform
from importlib.metadata import distribution, version
from PyInstaller.utils.hooks import collect_data_files, copy_metadata

if platform.system() != 'Darwin' or platform.machine() != 'arm64':
    raise SystemExit('QCBED.app must be built with native macOS arm64 Python')
base = Path(SPECPATH)
root = base.parent.parent
packages = ['qcbed','numpy','scipy','h5py','matplotlib','psutil','PyYAML',
            'PySide6','PySide6_Essentials','PySide6_Addons','shiboken6',
            'contourpy','cycler','fonttools','kiwisolver','packaging','pillow',
            'pyparsing','python-dateutil','six']
data = collect_data_files('qcbed')
# Include distribution notices/licenses and version metadata, not tests/tools.
for name in packages:
    data += copy_metadata(name)
    dist = distribution(name)
    for entry in dist.files or ():
        if any(word in str(entry).lower() for word in ('license','copying','notice')):
            file = Path(dist.locate_file(entry))
            if file.is_file(): data.append((str(file), 'licenses/'+name+'/'+str(entry.parent)))
a = Analysis([str(base/'launcher.py')],
    pathex=[str(root/'src'),str(base)], binaries=[], datas=data,
    hiddenimports=['scipy','scipy.linalg','scipy.fft','h5py','psutil','yaml','app_smoke'],
    hookspath=[], hooksconfig={'matplotlib': {'backends':['Agg', 'QtAgg']}},
    excludes=['pandas','sympy','torch','tensorflow','jax','cupy','numba','numexpr',
              'dask','xarray','sklearn','seaborn','pytest','_pytest','hypothesis','IPython','notebook','tkinter',
              'PyQt5','PyQt6','PySide2','numpy.tests','scipy.tests','h5py.tests','matplotlib.tests'],
    noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz,a.scripts,[],exclude_binaries=True,name='QCBED',
    debug=False,bootloader_ignore_signals=False,strip=False,upx=False,
    console=False,argv_emulation=False,target_arch='arm64',codesign_identity=None,entitlements_file=None)
collection = COLLECT(exe,a.binaries,a.datas,strip=False,upx=False,name='QCBED')
app = BUNDLE(collection,name='QCBED.app',bundle_identifier='org.qcbed.QCBED',
    info_plist={'CFBundleDisplayName':'QCBED','CFBundleShortVersionString':version('qcbed'),
                'CFBundleVersion':version('qcbed'),'NSHighResolutionCapable':True})

"""Platform reporting tests; never run benchmark simulations or timings."""
import builtins
from contextlib import nullcontext
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def benchmark(monkeypatch):
    # The standalone script imports a sibling and sets worker thread variables.
    # Restore those import-time environment changes after each test.
    directory = Path(__file__).resolve().parents[1] / "benchmarks"
    monkeypatch.syspath_prepend(str(directory))
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS"):
        monkeypatch.setenv(name, "1")
    spec = importlib.util.spec_from_file_location(
        "cpu_performance_audit", directory / "cpu_performance.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("system", ["darwin", "linux"])
def test_non_windows_does_not_import_registry(benchmark, monkeypatch, system):
    monkeypatch.setattr(benchmark.sys, "platform", system)
    monkeypatch.setattr(benchmark.platform, "processor", lambda: "")
    monkeypatch.setattr(benchmark.platform, "machine", lambda: "arm64")
    original = builtins.__import__

    def checked_import(name, *args, **kwargs):
        assert name != "winreg"
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", checked_import)
    assert benchmark._cpu_model() == "arm64"


def test_windows_registry_cpu_name(benchmark, monkeypatch):
    monkeypatch.setattr(benchmark.sys, "platform", "win32")
    registry = SimpleNamespace(
        HKEY_LOCAL_MACHINE=1, OpenKey=lambda *args: nullcontext("key"),
        QueryValueEx=lambda *args: ("  Example CPU  ", 1))
    monkeypatch.setitem(benchmark.sys.modules, "winreg", registry)
    assert benchmark._cpu_model() == "Example CPU"


@pytest.mark.parametrize("missing", [True, False])
def test_windows_registry_unavailable(benchmark, monkeypatch, missing):
    monkeypatch.setattr(benchmark.sys, "platform", "win32")
    monkeypatch.setattr(benchmark.platform, "processor", lambda: "fallback CPU")

    def denied(*args):
        raise PermissionError("registry denied")

    registry = None if missing else SimpleNamespace(HKEY_LOCAL_MACHINE=1, OpenKey=denied)
    monkeypatch.setitem(benchmark.sys.modules, "winreg", registry)
    assert benchmark._cpu_model() == "fallback CPU"


def test_library_paths_filtered_sorted_and_deduplicated(benchmark, monkeypatch):
    paths = ["/lib/libopenblas.so", "C:/libs/mkl_rt.dll", "/lib/unrelated.so",
             "/System/Accelerate.framework/Accelerate", "/lib/libopenblas.so"]
    process = SimpleNamespace(memory_maps=lambda: [SimpleNamespace(path=p) for p in paths])
    monkeypatch.setattr(benchmark.psutil, "Process", lambda: process)
    libraries, status = benchmark._loaded_math_libraries()
    assert libraries == sorted(set(paths) - {"/lib/unrelated.so"})
    assert status == "available"


@pytest.mark.parametrize("failure", ["missing", "unsupported", "denied", "os_error"])
def test_unavailable_memory_maps(benchmark, monkeypatch, failure):
    errors = {"unsupported": NotImplementedError("not supported"),
              "denied": benchmark.psutil.AccessDenied(pid=123),
              "os_error": OSError("inspection failed")}

    def mappings():
        raise errors[failure]

    process = SimpleNamespace() if failure == "missing" else SimpleNamespace(memory_maps=mappings)
    monkeypatch.setattr(benchmark.psutil, "Process", lambda: process)
    libraries, status = benchmark._loaded_math_libraries()
    assert libraries == []
    assert status.startswith("unavailable:")


def test_empty_supported_maps_are_not_reported_as_unavailable(benchmark, monkeypatch):
    monkeypatch.setattr(benchmark.psutil, "Process",
                        lambda: SimpleNamespace(memory_maps=lambda: []))
    assert benchmark._loaded_math_libraries() == ([], "available")


def test_environment_report_without_memory_maps(benchmark, monkeypatch):
    monkeypatch.setattr(benchmark.psutil, "Process", lambda: SimpleNamespace())
    report = benchmark.environment()
    assert report["cpu_model"]
    assert report["machine"] == benchmark.platform.machine()
    assert report["loaded_math_libraries"] == []
    assert report["loaded_math_libraries_status"].startswith("unavailable:")
    assert report["numpy_version"] == benchmark.np.__version__
    assert report["interpreter"] == benchmark.sys.executable
    assert "numpy_config" in report
    json.dumps(report)

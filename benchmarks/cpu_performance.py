"""Real-case profiling only. No production changes, installations or persistent pools.

Run with qcbed-python, from any directory. Results/checkpoints go to tmp/cpu_performance.
FFT substitutions are scoped to this benchmark process and preserve signs/norms.
Each measured configuration executes the complete existing production tilt engine.
"""
import os
INCOMING_ENV = {k: os.environ.get(k) for k in
                ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS')}
for _key in INCOMING_ENV:
    os.environ[_key] = '1'  # this process and spawned children only, before imports

from collections import defaultdict
from contextlib import contextmanager, redirect_stdout
from dataclasses import replace
import argparse
import cProfile
import importlib.util
import inspect
import io
import json
from pathlib import Path
import pickle
import platform
import pstats
import sys
import threading
from time import perf_counter, process_time
import tracemalloc
from unittest.mock import patch

import numpy as np
import psutil
import qcbed.refinement_trial as trial
import qcbed.tilt_series as series
import qcbed.parallel as parallel
from qcbed.scan import generate_scan_tilts
from qcbed.refinement import pack_refinement_parameters
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records
from qcbed.experimental import parse_experimental_intensities
from profile_thickness_only import thickness_only, profile_case

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'tmp/cpu_performance'
NP_FFT, NP_IFFT = np.fft.fft2, np.fft.ifft2


def residual(result, baseline):
    a, intensity = result.amplitudes, result.intensities
    ba, bi = baseline
    return dict(amplitude_max_abs=float(np.max(abs(a-ba))),
                intensity_max_abs=float(np.max(abs(intensity-bi))),
                intensity_relative_l2=float(np.linalg.norm(intensity-bi)/np.linalg.norm(bi)),
                bitwise_equal=bool(np.array_equal(a,ba) and np.array_equal(intensity,bi)))


def timed_run(fn):
    """100-ms RSS sampling; includes child working sets but not shared-page dedup."""
    proc = psutil.Process()
    start_rss = proc.memory_info().rss
    peak = [start_rss]; stop = threading.Event()
    def sample():
        while not stop.is_set():
            try:
                rss = proc.memory_info().rss
                for child in proc.children(recursive=True):
                    try: rss += child.memory_info().rss
                    except psutil.Error: pass
                peak[0] = max(peak[0],rss)
            except psutil.Error: pass
            stop.wait(.1)
    thread = threading.Thread(target=sample,daemon=True); thread.start()
    begin = perf_counter(); cpu = process_time()
    try: result = fn()
    finally:
        wall = perf_counter()-begin; cpu = process_time()-cpu
        stop.set(); thread.join()
    return result, dict(wall_seconds=wall,parent_cpu_seconds=cpu,
                       rss_start_mib=start_rss/2**20,peak_rss_mib=peak[0]/2**20,
                       peak_rss_increment_mib=(peak[0]-start_rss)/2**20)


@contextmanager
def fft_backend(name='numpy', workers=1):
    totals=defaultdict(float); counts=defaultdict(int)
    if name=='numpy': forward,inverse=NP_FFT,NP_IFFT
    elif name=='scipy':
        import scipy.fft
        forward,inverse=scipy.fft.fft2,scipy.fft.ifft2
    else: raise ValueError(name)
    def wrap(fn,label):
        def measured(*a,**kw):
            if name=='scipy':kw['workers']=workers
            t=perf_counter()
            try:return fn(*a,**kw)
            finally:totals[label]+=perf_counter()-t;counts[label]+=1
        return measured
    with patch.object(np.fft,'fft2',wrap(forward,'forward')), patch.object(np.fft,'ifft2',wrap(inverse,'inverse')):
        yield totals,counts


def _cpu_model():
    """Prefer the Windows registry when available; otherwise report architecture."""
    if sys.platform == 'win32':
        try:
            import winreg
            with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,r'HARDWARE\DESCRIPTION\System\CentralProcessor\0') as key:
                model=winreg.QueryValueEx(key,'ProcessorNameString')[0].strip()
            if model:
                return model
        except (ImportError, OSError):
            pass
    return platform.processor() or platform.machine() or 'unknown'


def _loaded_math_libraries():
    """Return detected paths and an explicit status for best-effort inspection."""
    try:
        mappings=psutil.Process().memory_maps()
    except (AttributeError, NotImplementedError, OSError, psutil.Error) as exc:
        return [], f'unavailable: {type(exc).__name__}: {exc}'
    paths=sorted({entry.path for entry in mappings
        if any(token in entry.path.lower() for token in
               ('mkl','blas','pocketfft','fftw','iomp','accelerate','veclib'))})
    return paths, 'available'


def environment():
    out=dict(interpreter=sys.executable,python=sys.version,numpy=np.__version__,
             platform=platform.platform(),logical_cpus=os.cpu_count(),physical_cpus=psutil.cpu_count(logical=False),
             inherited_environment=INCOMING_ENV,benchmark_environment={k:os.environ[k] for k in INCOMING_ENV},
             optional_packages={n:bool(importlib.util.find_spec(n)) for n in ('scipy','mkl_fft','pyfftw','threadpoolctl')})
    out['cpu_model']=_cpu_model()
    out['machine']=platform.machine()
    for name in ('numpy','scipy'):
        if importlib.util.find_spec(name):
            mod=__import__(name); buff=io.StringIO()
            with redirect_stdout(buff):mod.show_config()
            out[name+'_config']=buff.getvalue();out[name+'_version']=mod.__version__
    out['numpy_fft_module']=NP_FFT.__module__
    import numpy.fft._pocketfft as pocket
    out['numpy_fft_source']=pocket.__file__
    if importlib.util.find_spec('scipy'):
        import scipy.fft._pocketfft as spocket
        out['scipy_fft_source']=spocket.__file__
    # Trigger BLAS loading; this is not inside any measured simulation.
    np.ones((32,32)) @ np.ones((32,32))
    out['loaded_math_libraries'],out['loaded_math_libraries_status']=_loaded_math_libraries()
    out['psutil_version']=psutil.__version__
    for dll in out['loaded_math_libraries']:
        if Path(dll).name.lower().startswith('mkl_rt'):
            import ctypes
            try:
                lib=ctypes.CDLL(dll)
                fn=lib.MKL_Get_Max_Threads;fn.restype=ctypes.c_int;fn.argtypes=[]
                out['effective_mkl_max_threads']=fn()
            except (OSError,AttributeError):out['effective_mkl_max_threads']='query unavailable'
    out['conda_math_packages']=[]
    for f in (Path(sys.prefix)/'conda-meta').glob('*.json'):
        if any(t in f.name.lower() for t in ('mkl','blas','fft','numpy','scipy')):
            v=json.loads(f.read_text());out['conda_math_packages'].append({k:v.get(k) for k in ('name','version','build')})
    return out


def load_case():
    source=ROOT/'reference_data/FOR007.DAT'
    state,_=parse_legacy_refinement(source,seed=2026)
    state,_=thickness_only(state)
    records=_CaptureRecords(source);_parse_for007_records(records)
    row=records.values['output beam count and indices']
    beams=tuple(zip(row[1::2],row[2::2]))
    exp=parse_experimental_intensities(source.parent/'in.dat',beam_indices_hk=beams,scan=state.scan)
    return state,pack_refinement_parameters(state),exp


def baseline_case(state,vector,exp):
    evidence={}
    prepare,sequence=trial.prepare_multislice_sample,trial.simulate_slice_sequence
    def capture_prepare(*a,**kw):
        result=prepare(*a,**kw);evidence['prepared']=result;return result
    def capture_sequence(*a,**kw):
        if 'gratings' not in evidence:evidence['gratings']=a[0]
        return sequence(*a,**kw)
    with patch.object(trial,'prepare_multislice_sample',capture_prepare), patch.object(trial,'simulate_slice_sequence',capture_sequence):
        (result,profile),wall=timed_run(lambda:profile_case(state,vector,exp))
    prepared=evidence['prepared']
    # Ensure the isolated production sample is EXACTLY the trial's final grating.
    for stored,grating in zip(prepared.slice_types,evidence['gratings']):
        np.testing.assert_array_equal(stored.phase_grating,grating)
    return result,prepared,dict(**wall,profile=profile)


def line_profile(prepared,tilts,batch):
    """Time existing _chunk source lines, without replacing the propagation loop."""
    source,first=inspect.getsourcelines(series._chunk)
    line_text={first+i:s.strip() for i,s in enumerate(source)}
    loop_line=next(n for n,s in line_text.items() if s.startswith('for occurrence,'))
    totals=defaultdict(float); counts=defaultdict(int); previous={}
    def label(line):
        s=line_text.get(line,'')
        if line<loop_line:return 'kernel_setup_and_initialization'
        if 'wave = wave *' in s:return 'phase_grating_multiplication'
        if 'reciprocal *= kernels' in s:return 'propagator_multiplication'
        if 'batched_wave_real_to_reciprocal(' in s:return 'forward_dispatch_inclusive'
        if 'batched_wave_reciprocal_to_real(' in s:return 'inverse_dispatch_inclusive'
        if '_batch_array(' in s:return 'post_product_validation'
        if 'return reciprocal[' in s:return 'beam_extraction'
        return 'loop_control'
    def trace(frame,event,arg):
        if frame.f_code is not series._chunk.__code__:return None
        now=perf_counter();key=id(frame)
        if event in ('line','return'):
            if key in previous:
                line,t=previous[key];name=label(line);totals[name]+=now-t;counts[name]+=1
            if event=='return':previous.pop(key,None)
            else:previous[key]=(frame.f_lineno,perf_counter())
        return trace
    with fft_backend() as (fft,fft_counts):
        t=perf_counter();sys.settrace(trace)
        try:result=series.simulate_tilt_series_batched(prepared,tilts,batch_size=batch)
        finally:sys.settrace(None)
        wall=perf_counter()-t
    return result,dict(wall_seconds=wall,line_seconds=dict(totals),line_counts=dict(counts),
                       fft_seconds=dict(fft),fft_counts=dict(fft_counts),
                       unaccounted_seconds=wall-sum(totals.values()))


def _timed_initializer(prepared,batch):
    t=perf_counter();parallel._initialize_worker(prepared,batch)
    (DEST/f'worker_init_{os.getpid()}.json').write_text(json.dumps(dict(start=t,end=perf_counter(),pid=os.getpid())))


def _timed_chunk(start,stop,tilts):
    t=perf_counter();cpu=process_time()
    result=parallel._worker_chunk(start,stop,tilts)
    (DEST/f'worker_task_{os.getpid()}_{start}.json').write_text(json.dumps(dict(start=t,end=perf_counter(),cpu_seconds=process_time()-cpu,pid=os.getpid(),first=start,stop=stop)))
    return result


def process_profile(prepared,tilts,batch,workers):
    # New diagnostic files are isolated by deleting only this script's old markers.
    for pattern in ('worker_init_*.json','worker_task_*.json'):
        for p in DEST.glob(pattern):p.unlink()
    real=parallel.ProcessPoolExecutor
    timings=defaultdict(float)
    class TimedExecutor(real):
        def __init__(self,*a,**kw):
            kw['initializer']=_timed_initializer;t=perf_counter()
            super().__init__(*a,**kw);timings['executor_construction_seconds']+=perf_counter()-t
        def submit(self,fn,*a,**kw):
            t=perf_counter();f=super().submit(_timed_chunk,*a,**kw)
            timings['submission_seconds']+=perf_counter()-t
            original=f.result
            def get(*a,**kw):
                t=perf_counter()
                try:return original(*a,**kw)
                finally:timings['future_result_seconds']+=perf_counter()-t
            f.result=get
            return f
        def shutdown(self,*a,**kw):
            t=perf_counter()
            try:return super().shutdown(*a,**kw)
            finally:timings['shutdown_seconds']+=perf_counter()-t
    begin=perf_counter()
    with patch.object(parallel,'ProcessPoolExecutor',TimedExecutor):
        result,wall=timed_run(lambda:parallel.simulate_tilt_series_parallel(prepared,tilts,
            max_workers=workers,process_chunk_size=batch,batch_size=batch))
    tasks=[json.loads(p.read_text()) for p in DEST.glob('worker_task_*.json')]
    initializers=[json.loads(p.read_text()) for p in DEST.glob('worker_init_*.json')]
    timings['start_to_first_initializer_seconds']=min(x['start'] for x in initializers)-begin
    timings['initializer_body_sum_seconds']=sum(x['end']-x['start'] for x in initializers)
    timings['worker_computation_sum_seconds']=sum(x['end']-x['start'] for x in tasks)
    timings['worker_compute_span_seconds']=max(x['end'] for x in tasks)-min(x['start'] for x in tasks)
    timings['worker_cpu_sum_seconds']=sum(x['cpu_seconds'] for x in tasks)
    return result,dict(**wall,workers=workers,batch_size=batch,process_chunk_size=batch,
                       actual_workers=result.parallel.actual_workers,timings=dict(timings),tasks=tasks,initializers=initializers)


def microbenchmark(batch,backends):
    rng=np.random.default_rng(2026)
    a=(rng.normal(size=(batch,16,16))+1j*rng.normal(size=(batch,16,16))).astype(np.complex128)
    out=[]
    for name,workers in backends:
        if name=='numpy':functions=(NP_FFT,NP_IFFT)
        else:
            import scipy.fft
            functions=(scipy.fft.fft2,scipy.fft.ifft2)
        for direction,fn in zip(('forward','inverse'),functions):
            kwargs=dict(axes=(-2,-1),norm='forward')
            if name=='scipy':kwargs['workers']=workers
            expected=(NP_FFT if direction=='forward' else NP_IFFT)(a,axes=(-2,-1),norm='forward')
            actual=fn(a,**kwargs);np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=2e-12)
            reps=100;t=perf_counter()
            for _ in range(reps):fn(a,**kwargs)
            reps=max(100,int(.5/max((perf_counter()-t)/reps,1e-9)))
            values=[]
            for _ in range(3):
                t=perf_counter()
                for i in range(reps):fn(a,**kwargs)
                values.append((perf_counter()-t)/reps)
            per=float(np.median(values))
            out.append(dict(backend=name,workers=workers,direction=direction,shape=list(a.shape),
                            seconds_per_call=per,calls_per_second=1/per,repetitions=reps,
                            repeats=values,out_of_place=not np.shares_memory(a,actual),
                            max_abs_error=float(np.max(abs(actual-expected)))))
    return out


def scalar_diagnostics(prepared,tilts,baseline):
    """Small real-case prefix to locate scalar overhead, not a full-run estimate."""
    import qcbed.multislice as multislice
    tilts=tilts[:32];reference=(baseline[0][:32],baseline[1][:32])
    _,plain=timed_run(lambda:series.simulate_tilt_series_serial(prepared,tilts))
    times=defaultdict(float)
    def measured(fn,label):
        def wrapped(*a,**kw):
            t=perf_counter()
            try:return fn(*a,**kw)
            finally:times[label]+=perf_counter()-t
        return wrapped
    prof=cProfile.Profile()
    with patch.object(multislice,'apply_phase_grating',measured(multislice.apply_phase_grating,'phase_product_inclusive')), \
         patch.object(multislice,'apply_propagator',measured(multislice.apply_propagator,'propagator_product_inclusive')), \
         fft_backend() as (fft,counts):
        t=perf_counter();prof.enable()
        result=series.simulate_tilt_series_serial(prepared,tilts)
        prof.disable();wall=perf_counter()-t
    stream=io.StringIO();stats=pstats.Stats(prof,stream=stream)
    stats.strip_dirs().sort_stats('tottime').print_stats(40)
    (DEST/'scalar_prefix_profile.txt').write_text(stream.getvalue())
    selected=[]
    for (file,line,name),(cc,nc,tt,ct,callers) in stats.stats.items():
        if name in ('_array','_product','_centered_bins','reciprocal_index_grid',
                    'centered_to_unshifted','unshifted_to_centered','simulate_slice_sequence'):
            selected.append(dict(file=file,line=line,function=name,calls=nc,self_seconds=tt,inclusive_seconds=ct))
    return dict(tilts=len(tilts),plain=plain,profile_wall_seconds=wall,
                product_seconds=dict(times),fft_seconds=dict(fft),functions=selected,**residual(result,reference))


def layout_evidence(batch):
    a=np.ones((batch,16,16),dtype=np.complex128)
    import scipy.fft
    out={}
    for name,fn in [('numpy',NP_FFT),('scipy',scipy.fft.fft2)]:
        b=fn(a,axes=(-2,-1),norm='forward')
        out[name]=dict(shape=list(b.shape),input_strides=list(a.strides),output_strides=list(b.strides),
                       c_contiguous=bool(b.flags.c_contiguous),f_contiguous=bool(b.flags.f_contiguous))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume',action='store_true',help='reuse only completed stages in this output directory')
    args=parser.parse_args();DEST.mkdir(parents=True,exist_ok=True)
    report=json.loads((DEST/'results.json').read_text()) if args.resume and (DEST/'results.json').exists() else {}
    def save(label,value):
        report[label]=value
        (DEST/'results.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(label,json.dumps(value if label not in ('environment','microbenchmark') else {'saved':True}),flush=True)
    save('environment',environment())
    state,vector,exp=load_case();tilts=generate_scan_tilts(state.scan)
    if 'baseline' not in report:
        print('START full masked production objective',flush=True)
        result,prepared,metrics=baseline_case(state,vector,exp)
        np.savez(DEST/'baseline.npz',amplitudes=result.amplitudes,intensities=result.calculated_intensities)
        (DEST/'prepared.pkl').write_bytes(pickle.dumps(prepared,protocol=pickle.HIGHEST_PROTOCOL))
        save('baseline',metrics)
    prepared=pickle.loads((DEST/'prepared.pkl').read_bytes())
    with np.load(DEST/'baseline.npz') as data:baseline=(data['amplitudes'].copy(),data['intensities'].copy())
    save('configuration',dict(mesh=[prepared.mesh.mesh_y,prepared.mesh.mesh_x],tilts=len(tilts),
        beams=len(prepared.parameters.beam_indices_hk),occurrences=len(prepared.parameters.slice_sequence_zero_based),
        stored_slice_types=len(prepared.slice_types),depth_angstrom=float(vector[0]),objective_mode='masked',
        prepared_pickle_bytes=len(pickle.dumps(prepared,protocol=pickle.HIGHEST_PROTOCOL))))
    for batch in (1,4,8,16,32,64,128,256,512):
        key=f'batch_{batch}'
        if key in report:continue
        print('START',key,flush=True)
        result,metrics=timed_run(lambda:series.simulate_tilt_series_batched(prepared,tilts,batch_size=batch))
        metrics.update(batch_size=batch,ms_per_tilt=metrics['wall_seconds']*1000/len(tilts),**residual(result,baseline))
        np.testing.assert_allclose(result.amplitudes,baseline[0],rtol=2e-12,atol=2e-12)
        save(key,metrics)
    best=min((1,4,8,16,32,64,128,256,512),key=lambda b:report[f'batch_{b}']['wall_seconds'])
    save('best_batch_size',best)
    backends=[('numpy',1)]+([('scipy',w) for w in (1,2,4)] if importlib.util.find_spec('scipy') else [])
    for name,workers in backends:
        key=f'backend_{name}_{workers}'
        if key in report:continue
        print('START',key,flush=True)
        with fft_backend(name,workers) as (fft,counts):
            result,metrics=timed_run(lambda:series.simulate_tilt_series_batched(prepared,tilts,batch_size=best))
        metrics.update(fft_seconds=dict(fft),fft_counts=dict(counts),**residual(result,baseline))
        np.testing.assert_allclose(result.amplitudes,baseline[0],rtol=2e-12,atol=2e-12)
        save(key,metrics)
    if 'microbenchmark' not in report:save('microbenchmark',microbenchmark(best,backends))
    if 'line_profile' not in report:
        result,metrics=line_profile(prepared,tilts,best);metrics.update(residual(result,baseline));save('line_profile',metrics)
    if 'allocations' not in report:
        tracemalloc.start();t=perf_counter()
        result=series.simulate_tilt_series_batched(prepared,tilts,batch_size=best)
        wall=perf_counter()-t;current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
        save('allocations',dict(wall_seconds=wall,traced_current_bytes=current,traced_peak_bytes=peak,**residual(result,baseline)))
    if 'cprofile' not in report:
        prof=cProfile.Profile();t=perf_counter();prof.enable()
        result=series.simulate_tilt_series_batched(prepared,tilts,batch_size=best)
        prof.disable();wall=perf_counter()-t;prof.dump_stats(str(DEST/'batched.prof'))
        stream=io.StringIO();pstats.Stats(prof,stream=stream).strip_dirs().sort_stats('tottime').print_stats(40)
        (DEST/'batched_profile.txt').write_text(stream.getvalue())
        save('cprofile',dict(wall_seconds=wall,**residual(result,baseline)))
    for workers in (2,4):
        key=f'parallel_{workers}'
        if key in report:continue
        print('START',key,flush=True)
        result,metrics=process_profile(prepared,tilts,best,workers)
        metrics.update(residual(result,baseline));np.testing.assert_allclose(result.amplitudes,baseline[0],rtol=2e-12,atol=2e-12)
        save(key,metrics)
    if 'serialization' not in report:
        t=perf_counter()
        for _ in range(100):data=pickle.dumps(prepared,protocol=pickle.HIGHEST_PROTOCOL)
        dump=(perf_counter()-t)/100;t=perf_counter()
        for _ in range(100):pickle.loads(data)
        save('serialization',dict(bytes=len(data),dump_seconds=dump,load_seconds=(perf_counter()-t)/100))
    if 'serial_plain' not in report:
        print('START full minimally instrumented serial forward',flush=True)
        result,metrics=timed_run(lambda:series.simulate_tilt_series_serial(prepared,tilts))
        metrics.update(residual(result,baseline));save('serial_plain',metrics)
    for repeat in (1,2):
        key=f'best_repeat_{repeat}'
        if key not in report:
            result,metrics=timed_run(lambda:series.simulate_tilt_series_batched(prepared,tilts,batch_size=best))
            metrics.update(residual(result,baseline));save(key,metrics)
    if 'scalar_prefix_profile' not in report:
        save('scalar_prefix_profile',scalar_diagnostics(prepared,tilts,baseline))
    if 'layout_evidence' not in report and importlib.util.find_spec('scipy'):
        save('layout_evidence',layout_evidence(best))
    for workers in (2,4):
        key=f'parallel_plain_{workers}'
        if key not in report:
            result,metrics=timed_run(lambda:parallel.simulate_tilt_series_parallel(prepared,tilts,
                max_workers=workers,process_chunk_size=best,batch_size=best))
            metrics.update(residual(result,baseline));save(key,metrics)
    print('COMPLETE',DEST/'results.json',flush=True)


if __name__=='__main__':main()

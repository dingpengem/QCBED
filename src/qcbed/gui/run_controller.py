"""Spawned canonical execution without a Qt dependency or alternate run semantics."""
import multiprocessing
import os
from pathlib import Path
import time
import traceback
from contextlib import redirect_stdout, redirect_stderr

from qcbed import prepare_run_snapshot, run_snapshot


from .readiness import readiness


class _Stream:
    def __init__(self, connection): self.connection = connection
    def write(self, text):
        if text: self.connection.send(('log', text))
        return len(text)
    def flush(self): pass


def _worker(snapshot, resume, connection, cancel_event=None):
    # Set inside the isolated process; the editor's environment is untouched.
    os.environ['MPLBACKEND'] = 'Agg'
    started = time.monotonic()
    effective = None
    modern = snapshot.compiled.crystallography is not None and bool(snapshot.compiled.crystallography.site_groups)
    def observe(event):
        nonlocal effective
        connection.send(('progress', event))
        if not modern and effective is None and cancel_event is not None and cancel_event.is_set():
            effective = dict(message='Cancellation became effective; canonical retained-boundary finalization follows',
                boundary=event['classification'],
                evaluation=event['global_evaluation_index'], stage=event['stage'])
            connection.send(('cancellation_effective', effective))
            # The existing interruption path saves the last confirmed boundary
            # and finalizes it. Never inject into propagation or an atomic write.
            raise KeyboardInterrupt('GUI cooperative cancellation')
    try:
        connection.send(('started', {'pid': os.getpid()}))
        directory = Path(snapshot.compiled.configuration.run.output_directory)
        if cancel_event is not None and cancel_event.is_set():
            effective = dict(message='Cancellation became effective before canonical execution; no checkpoint created',
                             boundary='before_run', evaluation=None, stage=None)
            connection.send(('cancellation_effective', effective))
            connection.send(('result', dict(status='cancelled',
                summary=dict(status='cancelled', reason='Cancelled before starting canonical execution',
                             output_directory=str(directory), evaluations=0),
                cancellation=effective, elapsed_seconds=time.monotonic()-started,
                checkpoint=(directory/'checkpoint_latest.json').exists(),
                final_state=(directory/'final_parameters.json').exists())))
            return
        with redirect_stdout(_Stream(connection)), redirect_stderr(_Stream(connection)):
            options = dict(resume=resume,progress_hook=observe)
            if modern:options['cancel_requested']=cancel_event.is_set if cancel_event is not None else None
            result = run_snapshot(snapshot,**options)
        if modern and result.status=='cancelled':
            effective=dict(message='Cancellation effective before a new grouped evaluation; completed journal preserved',boundary='completed_evaluation')
            connection.send(('cancellation_effective',effective))
        directory = Path(snapshot.compiled.configuration.run.output_directory)
        if effective is None and cancel_event is not None and cancel_event.is_set():
            connection.send(('log', 'Cancellation arrived after the last safe progress boundary; canonical execution finished.'))
        connection.send(('result', dict(status='cancelled' if effective is not None and result.status != 'failed' else result.status,
            summary=result.summary, cancellation=effective,
            elapsed_seconds=time.monotonic()-started,
            checkpoint=(directory/'checkpoint_latest.json').exists(),
            final_state=(directory/'final_parameters.json').exists())))
    except BaseException:
        connection.send(('failure', traceback.format_exc()))
    finally:
        connection.close()


class RunController:
    """Poll from an event loop. No threads, YAML, forceful cancellation, or Qt."""
    def __init__(self):
        self.state = 'Idle'
        self.process = None
        self.connection = None
        self.result = None
        self.snapshot = None
        self.started = None
        self._terminal = False
        self.cancel_event = None
        self.cancel_pending = False

    @property
    def active(self): return self.process is not None

    def start(self, config, *, base_directory, resume=False, replace_outputs=False):
        if self.active: raise RuntimeError('A run is already active')
        issues = readiness(config, base_directory, resume=resume)
        errors = [i.message for i in issues if i.severity == 'error']
        if errors: raise ValueError('\n'.join(errors))
        modern=config.crystallography is not None
        if modern:
            from .output_paths import validate_output_directory,modern_execution_config
            validate_output_directory(config.run.output_directory,base_directory,create=True)
            config=modern_execution_config(config,base_directory,resume=resume)
        self.snapshot = prepare_run_snapshot(config, base_directory=base_directory)
        backups=()
        if not resume and not modern:
            from .output_paths import existing_outputs,backup_outputs
            if existing_outputs(config,base_directory):
                if not replace_outputs:raise ValueError('Confirm replacement of previous outputs before Fresh Run')
                backups=backup_outputs(config,base_directory)
        context = multiprocessing.get_context('spawn')
        receiver, sender = context.Pipe(duplex=False)
        cancel_event = context.Event()
        process = context.Process(target=_worker, args=(self.snapshot, resume, sender, cancel_event))
        try: process.start()
        except BaseException:
            for path,backup in reversed(backups):backup.rename(path)
            receiver.close()
            sender.close()
            raise
        sender.close()
        self.process, self.connection = process, receiver
        self.cancel_event = cancel_event
        self.cancel_pending = False
        self.result = None
        self.state = 'Running'
        self.started = time.monotonic()
        self._terminal = False

    def request_cancel(self):
        if not self.active or self._terminal or self.cancel_pending or self.cancel_event is None:
            return False
        self.cancel_pending = True
        self.cancel_event.set()
        return True

    def poll(self):
        events = []
        if not self.active: return events
        # Bound work per GUI tick even during verbose production runs.
        eof = False
        for _ in range(200):
            if not self.connection.poll(): break
            try: kind, payload = self.connection.recv()
            except EOFError:
                eof = True
                break
            events.append((kind, payload))
            if kind == 'result':
                self.result = payload
                self.state = {'completed':'Completed', 'failed':'Failed', 'cancelled':'Cancelled'}.get(payload['status'], 'Interrupted')
                self._terminal = True
            elif kind == 'failure':
                self.state = 'Failed'
                self._terminal = True
        if not self.process.is_alive() and (eof or not self.connection.poll()):
            self._finish(events)
        elif not self.process.is_alive() and self._terminal:
            self._finish(events)
        return events

    def _finish(self, events):
        self.process.join()
        if not self._terminal:
            self.state = 'Failed'
            events.append(('failure', f'Worker exited without a result (exit code {self.process.exitcode}).'))
        self.connection.close()
        self.process.close()
        self.process = self.connection = None

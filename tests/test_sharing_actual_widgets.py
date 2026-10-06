"""Real widget route: edit, Validate, Share dialogs, Validate, independent, Validate."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import pytest


@pytest.mark.parametrize('mode,n',[('real',1),('complex',2),('magnitude',1)])
def test_actual_widget_state_layers(tmp_path,mode,n):
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    from test_projected_elastic import config
    spec=importlib.util.spec_from_file_location('sharing_state_smoke',Path(__file__).parents[1]/'packaging/macos/sharing_state_smoke.py')
    harness=importlib.util.module_from_spec(spec);spec.loader.exec_module(harness)
    c=config(tmp_path,mode)
    c=replace(c,crystallography=replace(c.crystallography,projected_elastic=(),slice_boundaries_fractional=(0.,1.)))
    app=QApplication.instance() or QApplication([]);w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'p.yml');w.refresh();w.show()
    w.workflow.accept(0,c);assert w.validate_stage(1);assert w.validate_c_part('simulation')
    try:
        states=harness.exercise(w,app,tmp_path,mode)
        by={s['name']:s for s in states}
        for name in ('independent_validated','made_independent_validated','independent_revalidated'):
            s=by[name];assert s['vector_length']==1+2*n,s
            assert s['displayed']==f'{1+2*n} committed'
        shared=by['shared_validated'];assert shared['vector_length']==1+n,shared
        assert shared['bounds_length']==1+n and shared['displayed']==f'{1+n} committed'
        assert shared['physical']=='3 committed'
        assert all(r['group']=='Audit pair' for r in shared['committed'])
        for name in ('shared_before_validation','independent_before_validation'):
            assert by[name]['receipt']=='None' and not by[name]['run_enabled']
        assert f'{1+2*n} → {1+n}' in by['shared_before_validation']['message']
        assert f'{1+n} → {1+2*n}' in by['independent_before_validation']['message']
        assert shared['writeback_error'] is None
        assert len(shared['writeback'])==2
        assert shared['writeback'][0]['value'][0]==shared['writeback'][1]['value'][0]
        assert shared['receipt']!='None'
        assert len(by['selected_before_share']['selected'])==2
    finally:
        w.pending=False;w.document.saved_config=w.document.config;w.close()

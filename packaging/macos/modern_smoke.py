"""Opt-in frozen modern acceptance, driven through the production Qt window."""
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback
from dataclasses import replace


def main(manifest):
    if manifest.get("editor_controls_only"):return editor_controls_main(manifest)
    if manifest.get("live_abc_only"):return live_abc_main(manifest)
    if manifest.get("live_only"):return live_main(manifest)
    if manifest.get("compact_only"):return compact_main(manifest)
    if manifest.get("workflow_only"):return workflow_main(manifest)
    import numpy as np
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication
    from qcbed import load_yaml_config,save_yaml_config,compile_config_from_sources
    from qcbed.config_editing import Id
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    from qcbed.gui.readiness import readiness
    from qcbed.gui.results import load_results
    from qcbed.experimental_config import load_experimental_data
    output=Path(manifest['output'])
    if not output.is_absolute():raise ValueError('Absolute smoke output required')
    output.mkdir(parents=True,exist_ok=False)
    app=QApplication([]);window=ConfigWindow();window.show()
    report=dict(status='running',frozen=bool(getattr(sys,'frozen',False)),executable=sys.executable,
                cwd=os.getcwd(),architecture=platform.machine(),checks=[])
    stage='open';start=time.monotonic();cancel_sent=False;failed=False;initial=None
    def record(name):
        report['checks'].append(name)
        (output/'acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
    def set_config(c,name):
        c=replace(c,run=replace(c.run,output_directory=str(output/name),output_filename=str(output/name/'final_out.dat')))
        window.document=EditorDocument(c,output/(name+'.yml'));window.refresh();window.validate_document()
    def inspect():
        window.show_results();r=window.results_view.results
        assert r and r.history and r.summary['best'] and r.summary['last_evaluated']
        assert window.results_actions['Show Results'].isEnabled()
        assert len(r.parameters)==2 and r.comparison is not None
        a=r.comparison;assert np.isnan(a['experimental'][~a['mask']]).all()
        assert window.results_view.beam_details.rowCount()==9
        assert window.results_view.modern_history.rowCount()>0
        window.results_view.history_canvas.draw();window.results_view.physical_canvas.draw()
        window.results_view.pattern_canvas.draw()
        return r
    def tick():
        nonlocal stage,cancel_sent,failed,initial
        try:
            if failed:
                if not window.controller.active:app.exit(1)
                return
            if time.monotonic()-start>300:raise TimeoutError('Modern frozen smoke timeout')
            if stage=='open':
                import ase
                from ase.io import read,write
                report['ase']=ase.__version__
                window.open_document(manifest['reference_yaml']);window.validate_document()
                c=window.document.config
                assert c.crystallography and c.crystallography.site_groups[0].id=='Al_fcc'
                assert window.site_group_table.rowCount()==1
                assert 'strict_b5_compatible: True' in window.modern_summary.text()
                assert ('crystallography','cbed_basis','g1_hkl') in window.widgets
                record('modern_config_structure_projection_site_groups')
                if manifest.get('abc'):
                    assert [window.navigation.topLevelItem(i).text(0) for i in range(3)]==['A. Data Prep','B. Atomic Model','C. Main QCBED Parameters']
                    assert [a.get_title() for a in window.atomic_axes]==['Front — zone-axis projection','Top','Left']
                    v=window.data_panel
                    assert v.overview[0].shape==(256,256)
                    ax=v.overview_figure.axes[0]
                    assert len(ax.patches)==9 and len(ax.texts)==9
                    assert (ax.lines[0].get_xdata()[0],ax.lines[0].get_ydata()[0])==(127,126)
                    for beam in (0,8):
                        v.beams.setCurrentIndex(beam)
                        for axis,values in zip(v.axes,(v.data.intensities,np.sqrt(v.data.variances),v.data.weights)):
                            np.testing.assert_array_equal(axis.images[0].get_array(),values[:,beam].reshape(v.data.scan_rows,v.data.scan_columns))
                    assert window.save_document(output/'abc.yml')
                    if not manifest.get('reopen_only'):
                        original=window.data_panel.data
                        index=int(np.flatnonzero(original.weights[:,0])[0])
                        y,x=divmod(index,original.scan_columns)
                        def mask(**kw):
                            v=window.data_panel;v.beams.setCurrentIndex(0)
                            for spin,value in zip(v.box,(x,y,x+1,y+1)):spin.setValue(value)
                            v.mask(**kw)
                        mask()
                        assert window.data_panel.data.weights[index,0]==0
                        mask(restore=True)
                        np.testing.assert_array_equal(window.data_panel.data.weights,original.weights)
                        mask();mask(reset=True)
                        np.testing.assert_array_equal(window.data_panel.data.weights,original.weights)
                        mask()
                        np.testing.assert_array_equal(window.data_panel.data.intensities,original.intensities)
                        np.testing.assert_array_equal(window.data_panel.data.variances,original.variances)
                        assert window.save_document(overwrite=True)
                        snapshot=window.document.config
                        window.open_document(output/'abc.yml')
                        assert window.document.config==snapshot and window.data_panel.data.weights[index,0]==0
                        window.new_run_paths()
                        assert window.run_action.isEnabled() and not window.resume_action.isEnabled()
                    record('abc_linked_views_overview_mask_restore_reset_persistence_three_projections_new_paths')

                # Exercise dynamically imported CIF reader and spacegroup resources too.
                atoms=read(c.crystallography.structure.path)
                write(output/'probe.cif',atoms)
                cif=replace(c,crystallography=replace(c.crystallography,structure=replace(c.crystallography.structure,path=str(output/'probe.cif'),format='cif')))
                assert compile_config_from_sources(cif).crystallography
                record('ase_extxyz_and_cif_readers')
                key=('crystallography','site_groups',Id('Al_fcc'),'debye_waller_b')
                window.widgets[key].setText('0.87');assert window.apply_edits()
                for name in ('Conventional','Gradient'):
                    window.widgets[('experiment','gradient_option')].setCurrentText(name);assert window.apply_edits()
                assert window.save_document(output/'saved.yml')
                saved=window.document.config;window.open_document(output/'saved.yml');assert window.document.config==saved
                window.widgets[key].setText('0.86');assert window.save_document(overwrite=True)
                assert window.save_document(output/'save-as.yml')
                saved=window.document.config;window.open_document(output/'save-as.yml');assert window.document.config==saved
                record('B_edit_modes_save_save_as_reload')
                from qcbed.gui.compatibility import compatibility_report
                assert compatibility_report(saved,base_directory=output)['errors']
                record('modern_legacy_export_blocked')
                if manifest.get('reopen_only'):finish();return
                initial=load_yaml_config(manifest['run_yaml'])
                compiled=compile_config_from_sources(initial);exp=load_experimental_data(compiled,base_directory=output)
                assert np.count_nonzero(exp.weights)==7281 and np.array_equal(np.unique(exp.weights[exp.weights!=0]),[1.])
                record('gradient_7281_binary_support')
                set_config(initial,'gradient');assert window.run_action.isEnabled();window.start_run();stage='cancel'
            elif stage=='cancel':
                r=window.results_view.results
                if window.controller.active and r and r.summary['latest_valid_evaluation'] and not cancel_sent:
                    window.cancel_run();cancel_sent=True;assert not window.cancel_action.isEnabled()
                if window.controller.active or window.timer.isActive():return
                assert cancel_sent and window.controller.state=='Cancelled',window.run_log.toPlainText()
                r=inspect();assert r.summary['partial']
                first=next(v for v in r.history if v.get('objective') is not None)
                assert first['objective']==1.1896911365294378 and first['ynorm']==13447206.106343817
                record('spawn_cancel_after_completed_evaluation_results_accessible')
                window.grab().save(str(output/'cancelled.png'))
                changed=replace(window.document.config,crystallography=replace(initial.crystallography,depth_angstrom=841.))
                assert any(i.code=='resume_incompatible' for i in readiness(changed,output,resume=True))
                record('ownership_mismatch_rejected')
                assert window.resume_action.isEnabled();window.start_run(True);stage='resume'
            elif window.controller.active or window.timer.isActive():return
            elif stage=='resume':
                assert window.controller.state=='Completed',window.run_log.toPlainText()
                r=inspect();assert r.summary['retained_final']
                assert {v['stage'] for v in r.history}=={'ANNEAL','AMOEBA'}
                (output/'gradient-log.txt').write_text(window.run_log.toPlainText())
                report['gradient_summary']=window.controller.result['summary']
                record('gradient_resume_anneal_amoeba_history_masked_patterns')
                c=load_yaml_config(manifest['conventional_yaml']);set_config(c,'conventional')
                window.start_run();stage='conventional'
            elif stage=='conventional':
                assert window.controller.state=='Completed',window.run_log.toPlainText()
                r=inspect();assert r.summary['signal_mode']=='conventional'
                assert int(r.comparison['mask'].sum())==7974
                record('conventional_run_results_masked_patterns')
                finish()
        except BaseException:
            report.update(status='failed',traceback=traceback.format_exc());record('failure')
            failed=True;window.cancel_run()
            if not window.controller.active:app.exit(1)
    def finish():
        report.update(status='passed',elapsed_seconds=time.monotonic()-start)
        timer.stop()
        app.processEvents()
        window.data_panel.canvas.draw();window.data_panel.overview_canvas.draw()
        window.grab().save(str(output/'window.png'));record('modern_frozen_complete')
        window.pending=False;window.document.saved_config=window.document.config
        timer.stop();window.close();app.quit()
    timer=QTimer();timer.setInterval(50);timer.timeout.connect(tick);timer.start()
    return app.exec()


def workflow_main(manifest):
    """GUI-only acceptance. Never start a runner or numerical worker."""
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import QTimer
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    import numpy as np
    output=Path(manifest['output']);output.mkdir(parents=True,exist_ok=False)
    app=QApplication([]);window=ConfigWindow(manifest['reference_yaml']);window.show()
    report=dict(status='running',frozen=bool(getattr(sys,'frozen',False)),checks=[])
    def check():
        try:
            assert window.workflow_tabs.count()==3
            assert window.save_document(output/'working.yml')
            c=window.document.config
            existing=output/'run';existing.mkdir();(existing/'old.txt').write_text('existing output')
            window.document=EditorDocument(replace(c,run=replace(c.run,output_directory=str(existing),output_filename=str(existing/'final.dat'))),output/'working.yml')
            window.refresh()
            v=window.data_panel;assert v.overview[0].shape==(256,256)
            old=v.data;v.beams.setCurrentIndex(0)
            v.select_feature.setChecked(True)
            from types import SimpleNamespace
            v.pick_overview(SimpleNamespace(xdata=127.,ydata=126.,inaxes=v.overview_figure.axes[0]))
            ax=window.data_panel.overview_figure.axes[0]
            assert len(ax.patches)==9 and len(ax.texts)==9
            assert ax.lines[0].get_label()=='Zone axis'
            index=int(np.flatnonzero(old.weights[:,0])[0]);y,x=divmod(index,old.scan_columns)
            v=window.data_panel
            for spin,value in zip(v.box,(x,y,x+1,y+1)):spin.setValue(value)
            v.mask();assert window.data_panel.data.weights[index,0]==0
            np.testing.assert_array_equal(old.intensities,window.data_panel.data.intensities)
            np.testing.assert_array_equal(old.variances,window.data_panel.data.variances)
            assert not window.run_action.isEnabled()
            for stage in range(3):assert window.validate_stage(stage),window.stage_status[stage].text()
            assert window.run_action.isEnabled() and not window.resume_action.isEnabled()
            assert existing.exists() and (existing/'old.txt').read_text()=='existing output'
            saved=window.document.config
            assert window.save_document(output/'saved.yml')
            window.open_document(output/'saved.yml');assert window.document.config==saved
            assert not window.run_action.isEnabled()
            for stage in range(3):assert window.validate_stage(stage)
            app.processEvents()
            window.data_panel.overview_figure.savefig(output/'overview.png')
            window.data_panel.figure.savefig(output/'linked.png')
            window.grab().save(str(output/'window.png'))
            report.update(status='passed',checks=['three_tabs','manual_feature_origin_marker','nine_beam_boxes_labels',
                'mask_only_weights','ordered_validation','existing_output_ready_without_mutation','save_reload'],refinement_executed=False)
            code=0
        except BaseException:
            report.update(status='failed',traceback=traceback.format_exc());code=1
        (output/'acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
        window.pending=False;window.document.saved_config=window.document.config;window.close();app.exit(code)
    QTimer.singleShot(100,check)
    return app.exec()

def initialize_only(config,base):
    """Exercise canonical output/checkpoint setup, stop before optimizer entry."""
    from unittest.mock import patch
    from qcbed import prepare_run_snapshot,run_snapshot
    from qcbed.grouped_run import load_grouped_checkpoint,grouped_history
    class SmokeBoundary(Exception):pass
    snapshot=prepare_run_snapshot(config,base_directory=base)
    with patch('qcbed.grouped_run._run_stage',side_effect=SmokeBoundary('intentional packaging initialization stop')) as stage:
        try:run_snapshot(snapshot)
        except SmokeBoundary:pass
        else:raise AssertionError('Initialization did not reach the expected stage boundary')
        assert stage.call_count==1
    out=Path(config.run.output_directory)
    cp=load_grouped_checkpoint(out,snapshot)
    assert grouped_history(out,cp)==[] and cp.get('evaluations',0)==0
    assert (out/'checkpoint_latest.json').is_file() and (out/'patterns').is_dir()
    assert not (out/'.run.lock').exists()
    return dict(evaluations=0,optimizer_entered=False,checkpoint_status=cp['status'],
                stop='intentional smoke boundary before optimizer',directory=str(out))


def compact_main(manifest):
    """Frozen compact frontend acceptance; no optimization or forward evaluation."""
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import QTimer
    from qcbed import load_yaml_config
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    from qcbed.gui.data_model import provenance
    import qcbed.modern_yaml
    import yaml
    import numpy as np
    output=Path(manifest['output']);output.mkdir(parents=True,exist_ok=False)
    app=QApplication([]);window=ConfigWindow();window.show()
    report=dict(status='running',frozen=bool(getattr(sys,'frozen',False)),cwd=os.getcwd(),
                architecture=platform.machine(),module_file=qcbed.modern_yaml.__file__,checks=[])
    def check():
        try:
            if report['frozen']:
                assert Path(qcbed.modern_yaml.__file__).is_relative_to(Path(sys._MEIPASS))
            window.open_document(manifest['reference_yaml'])
            reference=window.document.config
            assert window.workflow_tabs.count()==3 and window.site_group_table.rowCount()>0
            assert len(window.atomic_axes)==3 and window.refinement_table.rowCount()>0
            assert window.data_panel.data is not None
            report['checks'].append('canonical_fcc_al_compact_load')
            window.open_document(manifest['interactive_yaml'])
            c=window.document.config;meta=provenance(c);data=window.data_panel.data
            assert meta['excluded_indices'] and meta['zone_axis_selection']
            assert window.data_panel.overview[0].shape==(256,256)
            ax=window.data_panel.overview_figure.axes[0]
            assert len(ax.patches)==9 and len(ax.texts)==9 and ax.lines[0].get_label()=='Zone axis'
            origin=meta['zone_axis_selection']['zone_axis_origin_px']
            assert [ax.lines[0].get_xdata()[0],ax.lines[0].get_ydata()[0]]==origin
            for tilt,beam in meta['excluded_indices']:assert data.weights[tilt,beam]==0
            report['checks'].append('external_arrays_pattern_mask_origin_provenance_annotations')
            assert window.save_document(output/'saved.yml')
            assert window.save_document(overwrite=True)
            assert window.save_document(output/'save-as.yml')
            saved=window.document.config
            for name in ('saved.yml','save-as.yml'):
                assert set(yaml.safe_load((output/name).read_text()))=={'schema_version','data_prep','atomic_model','qcbed'}
                window.open_document(output/name)
                assert window.document.config==saved
                assert provenance(saved)==meta
                for attr in ('intensities','variances','weights'):
                    np.testing.assert_array_equal(getattr(window.data_panel.data,attr),getattr(data,attr))
            report['checks'].append('compact_save_save_as_reload_exact')
            for path,modern in ((manifest['verbose_yaml'],True),(manifest['legacy_yaml'],False)):
                window.open_document(path)
                assert window.document.config==load_yaml_config(path)
                assert (window.document.config.crystallography is not None)==modern
                if modern:assert window.data_panel.data is not None
            report['checks'].append('verbose_modern_and_legacy_open')
            window.document=EditorDocument(replace(saved,run=replace(saved.run,
                output_directory=str(output/'initialization'),output_filename=str(output/'initialization/final_out.dat'))),output/'run.yml')
            window.refresh()
            for stage in range(3):assert window.validate_stage(stage),window.stage_status[stage].text()
            assert window.run_action.isEnabled()
            report['initialization']=initialize_only(window.document.config,output)
            report['checks'].append('canonical_run_checkpoint_output_initialization_zero_evaluations')
            app.processEvents();window.grab().save(str(output/'window.png'))
            report.update(status='passed',refinement_executed=False)
            code=0
        except BaseException:
            report.update(status='failed',traceback=traceback.format_exc());code=1
        (output/'acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
        window.pending=False
        if window.document:window.document.saved_config=window.document.config
        window.close();app.exit(code)
    QTimer.singleShot(100,check)
    return app.exec()


def live_main(manifest):
    """Spawn real Patch A from synthetic raw input; never enter refinement."""
    import numpy as np
    from types import SimpleNamespace
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import QTimer
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.data_model import read_data
    out=Path(manifest['output']);out.mkdir(parents=True,exist_ok=False)
    y,x=np.indices((180,180));raw=np.full((180,180),10.)
    for cx,cy in ((90,90),(140,90),(90,140)):
        raw+=100/(1+np.exp(np.clip((np.hypot(x-cx,y-cy)-20)*3,-700,700)))
    raw+=30*np.exp(-((x-90)**2+(y-90)**2)/10)
    pattern=out/'raw.npy';np.save(pattern,raw)
    aperture=100.*(((x-90)/20)**2+((y-90)/15)**2<1);np.save(out/'aperture.npy',aperture)
    app=QApplication([]);window=ConfigWindow();window.show();v=window.preparation
    report=dict(status='running',frozen=bool(getattr(sys,'frozen',False)),checks=[],refinement_executed=False)
    started=time.monotonic();launched=False
    def drag(a,b):
        v.begin(SimpleNamespace(xdata=a[0],ydata=a[1],inaxes=True))
        v.end(SimpleNamespace(xdata=b[0],ydata=b[1],inaxes=True))
    def finish(code):
        report['elapsed']=time.monotonic()-started
        (out/'acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
        timer.stop();window.pending=False
        if window.document:window.document.saved_config=window.document.config
        window.close();app.exit(code)
    def tick():
        nonlocal launched
        try:
            if not launched:
                assert v.mode.currentText()=='New Data Preparation'
                v.raw_path=out/'project.yml';v.inputs['pattern'].setText(str(pattern));v.load_pattern()
                v.inputs['calibration_aperture'].setText(str(out/'aperture.npy'));v.load_profile();v.accept_preview()
                assert v.geometry.radius_x_px==19.5 and v.geometry.radius_y_px==14.5
                # Aperture input/threshold UI checked; actual preparation uses
                # explicitly selected no-aperture path, not a pretend calibration.
                v.inputs['calibration_aperture'].clear()
                drag((70,90),(110,90));assert v.preview;v.accept_preview()
                v.tool.setCurrentIndex(1);drag((86.8,87.2),(92.7,93.1));v.accept_preview()
                v.binning.setValue(1)
                v.start();assert v.active and v.process.pid!=os.getpid()
                report['checks']+=['new_default','raw_input','optional_aperture_threshold_ellipse','manual_diameter','explicit_roi_snap','spawned_worker']
                launched=True;return
            if v.active:
                if time.monotonic()-started>120:raise RuntimeError('Preparation acceptance timeout')
                return
            assert window.document is not None,v.status.text()
            assert window.data_panel.data.intensities.shape==(1600,9)
            assert len(window.data_panel.overview_figure.axes[0].patches)==9
            assert not window.validate_stage(0)
            window.data_panel.select_zone();assert window.validate_stage(0),window.stage_status[0].text()
            before=read_data(window.document.config,window.document.base_directory)
            window.data_panel.mask()
            after=read_data(window.document.config,window.document.base_directory)
            np.testing.assert_array_equal(before.intensities,after.intensities)
            np.testing.assert_array_equal(before.variances,after.variances)
            assert window.validate_stage(0)
            assert v.mode.currentText()=='Load Previous Preparation'
            window.grab().save(str(out/'window.png'))
            report['checks']+=['live_beams','a8_refined_origin_confirmation','pattern_sigma_weight','weight_only_mask','validate_A_without_B','replay_mode']
            report['status']='passed';finish(0)
        except BaseException:
            report.update(status='failed',traceback=traceback.format_exc());finish(1)
    timer=QTimer();timer.setInterval(100);timer.timeout.connect(tick);timer.start()
    return app.exec()


def live_abc_main(manifest):
    """Real raw preparation through A/B/C; deliberate stop before optimizer entry."""
    import numpy as np
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    from qcbed.gui.preparation_model import prepare_files,new_preparation_document,apply_preparation
    from qcbed.gui.data_model import select_zone_axis,edit_mask,read_data,restored_overview
    from qcbed.data_preparation.live import LiveSelection
    from qcbed.data_preparation.models import DiskGeometry,SignalMode
    from qcbed.configuration import (StructureFileConfig,CbedReflectionBasisConfig,AtomicSiteGroup,
                                    ParameterBinding,TargetBound,BoundsConfig)
    from qcbed.execution import ForwardExecutionConfig
    from qcbed.crystallography.site_groups import group_target,requested_depth_target
    from qcbed import compile_config_from_sources
    out=Path(manifest['output']);out.mkdir(parents=True,exist_ok=False)
    app=QApplication.instance() or QApplication([]);window=ConfigWindow()
    report=dict(status='running',frozen=bool(getattr(sys,'frozen',False)),refinement_executed=False)
    try:
        geometry=DiskGeometry(*manifest['radius_px'])
        selection=LiveSelection(geometry,tuple(manifest['symmetry_roi']),76,76,final_binning=2,
            signal_mode=SignalMode.GRADIENT,symmetry_crop_half_size_px=(19,19))
        if manifest.get('interactive'):
            from types import SimpleNamespace
            v=window.preparation;window.show();v.raw_path=out/'project.yml'
            assert v.mode.currentText()=='New Data Preparation'
            assert [window.workflow_tabs.tabText(i) for i in range(3)]==['A. Data Prep','B. Atomic Model','C. Main QCBED Parameters']
            v.inputs['pattern'].setText(manifest['pattern']);v.load_pattern()
            def drag(a,b):
                v.begin(SimpleNamespace(xdata=a[0],ydata=a[1],inaxes=True))
                v.end(SimpleNamespace(xdata=b[0],ydata=b[1],inaxes=True))
            drag(*manifest['diameter_endpoints']);assert v.geometry is None
            v.accept_preview();assert v.geometry==geometry
            v.tool.setCurrentIndex(1)
            roi=manifest['symmetry_roi'];drag(roi[:2],roi[2:]);assert v.roi is None
            v.accept_preview();assert v.roi==tuple(roi)
            v.binning.setValue(2)
            v.signal.setCurrentText('gradient');v.run.click()
            assert v.active and v.process.pid!=os.getpid(),v.status.text()
            deadline=time.monotonic()+120
            while v.active and time.monotonic()<deadline:
                app.processEvents();time.sleep(.02)
            assert not v.active and window.document is not None,v.status.text()
            panel=window.data_panel;panel.select_zone();panel=window.data_panel
            assert len(panel.overview_figure.axes[0].patches)==9
            assert len(panel.overview_figure.axes[0].texts)==9
            assert any(line.get_label()=='Zone axis' for line in panel.overview_figure.axes[0].lines)
            for beam in range(9):
                panel.beams.setCurrentIndex(beam)
                assert [ax.get_title() for ax in panel.axes]==['Pattern','Sigma','Weight']
                expected=(panel.data.intensities,np.sqrt(panel.data.variances),panel.data.weights)
                for ax,array in zip(panel.axes,expected):
                    np.testing.assert_array_equal(np.asarray(ax.images[0].get_array()).ravel(),array[:,beam])
            panel.beams.setCurrentIndex(0)
            before=read_data(window.document.config,out)
            for spin,value in zip(panel.box,(18,18,19,19)):spin.setValue(value)
            panel.mask();config=window.document.config
            after=read_data(config,out)
            np.testing.assert_array_equal(before.intensities,after.intensities)
            np.testing.assert_array_equal(before.variances,after.variances)
            assert before.weights[18*38+18,0]==1 and after.weights[18*38+18,0]==0
            assert np.count_nonzero(before.weights!=after.weights)==1
            panel=window.data_panel;panel.overview_figure.savefig(out/'live-overview.png');panel.figure.savefig(out/'linked-pattern-sigma-weight.png')
            _,m=restored_overview(config,out)
            assert m['live_preparation']['correction']=='none'
            assert m['live_preparation']['noise_source']=='Poisson / intensity fallback'
            report['interactive_checks']=['manual_diameter_accept','centrosymmetric_roi_accept','spawned_live_worker',
                'nine_boxes_labels','zone_axis_marker','nine_linked_pattern_sigma_weight','one_pixel_weight_only']
        else:
            descriptor=prepare_files(manifest['pattern'],selection,out/'live')
            config=new_preparation_document(descriptor)
            config=select_zone_axis(config,out,beam=0)
            config=edit_mask(config,out,beam=0,box=(18,18,19,19))
            window.document=EditorDocument(config,out/'project.yml');window.refresh();window.show()
        assert window.validate_stage(0),window.stage_status[0].text()
        validated=window.workflow.preparation_descriptor
        replay=apply_preparation(config,validated)
        assert replay==config
        window.preparation.install(str(validated))
        assert window.document.config==config and window.preparation.mode.currentIndex()==1
        for name in ('intensities','variances','weights'):
            np.testing.assert_array_equal(getattr(read_data(replay,out),name),getattr(read_data(config,out),name))
        depth,b=requested_depth_target(),group_target('Al_fcc')
        config=replace(config,execution=ForwardExecutionConfig('batched',64),
            crystallography=replace(config.crystallography,
                structure=StructureFileConfig(manifest['structure'],'extxyz'),
                cbed_basis=CbedReflectionBasisConfig((2,0,0),(0,2,0)),depth_angstrom=840.,
                slice_boundaries_fractional=(0.,.5,1.),
                site_groups=(AtomicSiteGroup('Al_fcc','Al',(0,1,2,3),.86,1.,(0.,1.,1.)),)),
            refinement=replace(config.refinement,
                parameters=(ParameterBinding(depth,(1.05,.95),.5),ParameterBinding(b,(1.05,.95),.0002)),
                bounds=BoundsConfig('reject',(TargetBound(depth,600.,1200.),TargetBound(b,.5,1.5)))))
        window.document.commit_structure(config);window.refresh()
        assert window.validate_stage(1),window.stage_status[1].text()
        assert window.validate_stage(2),window.stage_status[2].text()
        assert window.run_action.isEnabled()
        assert [ax.get_title().split(' —')[0] for ax in window.atomic_axes]==['Front','Top','Left']
        assert window.site_group_table.rowCount()==1
        window.atomic_axes[0].figure.savefig(out/'atomic-views.png')
        compiled=compile_config_from_sources(config)
        mapping=[dict(lattice=list(b.experimental_index),hkl=list(b.original_hkl)) for b in compiled.crystallography.indexed_beams]
        assert all(row['hkl']==[2*row['lattice'][0],2*row['lattice'][1],0] for row in mapping)
        assert window.save_document(out/'saved.yml')
        saved=window.document.config;window.open_document(out/'saved.yml');assert window.document.config==saved
        _,metadata=restored_overview(config,out)
        if manifest.get('interactive'):
            import multiprocessing
            ctx=multiprocessing.get_context('spawn');receiver,sender=ctx.Pipe(duplex=False)
            process=ctx.Process(target=initialization_smoke_worker,args=(config,out,sender))
            process.start();sender.close();deadline=time.monotonic()+120
            while not receiver.poll() and process.is_alive() and time.monotonic()<deadline:
                app.processEvents();time.sleep(.02)
            assert receiver.poll(),'Initialization worker did not return'
            initialization=receiver.recv();receiver.close();process.join(10)
            assert process.exitcode==0 and 'error' not in initialization,initialization
            assert initialization['pid']!=os.getpid()
        else:initialization=initialize_only(config,out)
        report.update(initialization=initialization,
            corners={k:list(getattr(config.scan,k)) for k in ('corner_a','corner_b','corner_c')},
            reciprocal_frame=metadata['reciprocal_frame'],mapping=mapping,prepared_descriptor=str(validated),
            replay='exact canonical and array equality',validation='A/B/C passed',status='passed')
        window.grab().save(str(out/'window.png'))
        code=0
    except BaseException:
        report.update(status='failed',traceback=traceback.format_exc());code=1
    finally:
        (out/'acceptance.json').write_text(json.dumps(report,indent=2)+'\n')
        if window.document:window.document.saved_config=window.document.config
        window.pending=False;window.close();app.processEvents()
    return code


def initialization_smoke_worker(config,base,connection):
    """Test-only spawned canonical initialization, no optimizer execution."""
    try:connection.send(dict(initialize_only(config,base),pid=os.getpid()))
    except BaseException:connection.send(dict(error=traceback.format_exc()))
    finally:connection.close()


def editor_controls_main(manifest):
    """Focused frozen widget visibility check; no worker or objective execution."""
    from PySide6.QtWidgets import QApplication, QFileDialog
    from qcbed.gui.window import ConfigWindow
    from qcbed import load_yaml_config
    app=QApplication([])
    output=Path(manifest['output']);output.mkdir(parents=True,exist_ok=False)
    w=ConfigWindow(manifest['reference_yaml']);w.show();app.processEvents()
    try:
        QFileDialog.getOpenFileName=lambda *a,**k:(manifest['cif'],'')
        w.browse(w.widgets[('crystallography','structure','path')],'file')
        app.processEvents()
        table=w.site_group_table
        assert table.isVisible() and table.visibleRegion().boundingRect().height()>0
        assert table.item(0,0).text()=='Al_fcc'
        assert table.item(0,3).text()==''
        table.item(0,3).setText('0.86');w.apply_site_groups.click()
        app.processEvents()
        if manifest.get('refinement_drafts'):
            from PySide6.QtCore import Qt
            w.workflow_tabs.setCurrentIndex(2)
            w.navigation.setCurrentItem(w.section_items['Refinement'])
            w.widgets[('crystallography','depth_angstrom')].setText('1000')
            t=w.refinement_table
            rows={t.item(r,0).data(Qt.ItemDataRole.UserRole).category.value:r for r in range(t.rowCount())}
            t.item(rows['depth'],1).setCheckState(Qt.CheckState.Checked)
            assert not w.apply_edits()
            assert any('requires D1, D2, PLAMDA' in issue.message for issue in w.document.issues)
            for category in ('depth','debye_waller'):
                r=rows[category];t.item(r,1).setCheckState(Qt.CheckState.Checked)
                for col,value in ((4,'.8'),(5,'1.2'),(6,'.1'),
                        (8,'600' if category=='depth' else '.5'),
                        (9,'1200' if category=='depth' else '1.5')):
                    t.item(r,col).setText(value)
            assert w.apply_edits()
            assert w.document.config.crystallography.depth_angstrom==1000.
            if manifest.get('dynamic_refinement'):
                from qcbed import compile_config_from_sources
                from qcbed.grouped_run import check_grouped
                t=w.refinement_table
                for r in range(t.rowCount()):
                    target=t.item(r,0).data(Qt.ItemDataRole.UserRole)
                    if target.category.value=='depth':
                        t.item(r,1).setCheckState(Qt.CheckState.Unchecked)
                    if target.category.value=='scan_coordinate' and target.identifiers==('A','first'):
                        assert t.item(r,1).flags() & Qt.ItemFlag.ItemIsEnabled
                        t.item(r,1).setCheckState(Qt.CheckState.Checked)
                        for col,value in ((4,'.9'),(5,'1.1'),(6,'.01'),(8,'-2'),(9,'2')):
                            t.item(r,col).setText(value)
                assert w.apply_edits()
                assert [p['category'] for p in check_grouped(compile_config_from_sources(w.document.config))]==['debye_waller','scan_coordinate']
            app.processEvents()
            w.grab().save(str(output/'refinement-controls.png'))
        # Existing A input is a fixture; no preparation or numerical evaluation.
        w.workflow.accept(0,w.document.config);w.update_stages()
        w.stage_buttons[1].click()
        assert w.workflow.receipts[1] is not None,w.stage_status[1].text()
        w.workflow_tabs.setCurrentIndex(2)
        w.navigation.setCurrentItem(w.section_items['Numerical / Simulation cell'])
        app.processEvents()
        field=w.widgets[('crystallography','slice_boundaries_fractional')]
        assert field.isVisible() and field.visibleRegion().boundingRect().width()>100
        field.setText('0, 0.5, 1')
        w.stage_buttons[2].click()
        assert w.document.config.crystallography.slice_boundaries_fractional==(0.,.5,1.)
        if manifest.get('refinement_drafts'):
            assert w.document.files_checked
            assert 'Not checked' not in w.summary.toPlainText()
        assert w.save_document(output/'saved.yml')
        w.open_document(output/'saved.yml')
        w.workflow_tabs.setCurrentIndex(2)
        w.navigation.setCurrentItem(w.section_items['Numerical / Simulation cell'])
        app.processEvents()
        field=w.widgets[('crystallography','slice_boundaries_fractional')]
        assert field.isVisible()
        assert load_yaml_config(output/'saved.yml')==w.document.config
        w.grab().save(str(output/'visible-controls.png'))
        (output/'result.json').write_text(json.dumps(dict(status='PASS',frozen=bool(getattr(sys,'frozen',False)),
            slice_text=field.text(),groups=len(w.document.config.crystallography.site_groups),optimizer_evaluations=0)))
        return 0
    except Exception:
        (output/'failure.txt').write_text(traceback.format_exc());return 1
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()

"""Real raw-data handoff, never importing a historical preparation or running optimization."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pytest
from qcbed.data_preparation.live import LiveSelection
from qcbed.data_preparation.models import DiskGeometry,SignalMode
from qcbed.gui.preparation_model import prepare_files,new_preparation_document,apply_preparation
from qcbed.gui.data_model import select_zone_axis,edit_mask,read_data,restored_overview
from qcbed.gui.workflow import Workflow,validate_data
from qcbed.configuration import StructureFileConfig,CbedReflectionBasisConfig,AtomicSiteGroup
from qcbed import compile_config_from_sources,prepare_run_snapshot,run_snapshot

ROOT=Path(__file__).resolve().parents[1]


def atomic_model(config):
    # Only explicitly supplied B/C physical controls; no reference YAML/manifest.
    from qcbed.configuration import ParameterBinding,TargetBound,BoundsConfig
    from qcbed.crystallography.site_groups import group_target,requested_depth_target
    depth,b=requested_depth_target(),group_target('Al_fcc')
    from qcbed.execution import ForwardExecutionConfig
    return replace(config,execution=ForwardExecutionConfig('batched',64),refinement=replace(config.refinement,parameters=(ParameterBinding(depth,(1.05,.95),.5),ParameterBinding(b,(1.05,.95),.0002)),bounds=BoundsConfig('reject',(TargetBound(depth,600.,1200.),TargetBound(b,.5,1.5)))),crystallography=replace(config.crystallography,
        structure=StructureFileConfig(str(ROOT/'examples/fcc_al_conventional.extxyz'),'extxyz'),
        cbed_basis=CbedReflectionBasisConfig((2,0,0),(0,2,0)),
        depth_angstrom=840.,slice_boundaries_fractional=(0.,.5,1.),
        site_groups=(AtomicSiteGroup('Al_fcc','Al',(0,1,2,3),.86,1.,(0.,1.,1.)),)))


@pytest.fixture(scope='module')
def live_config(tmp_path_factory):
    base=tmp_path_factory.mktemp('live_abc')
    path=prepare_files(ROOT/'reference_data/pattern.dat',
        LiveSelection(DiskGeometry(38,38),(111,105,155,151),76,76,final_binning=2,
                      signal_mode=SignalMode.GRADIENT,symmetry_crop_half_size_px=(19,19)),base/'prep')
    c=select_zone_axis(new_preparation_document(path),base,beam=0)
    c=edit_mask(c,base,beam=0,box=(18,18,19,19))
    return c,base


def test_live_a_snapshot_replay_and_indexing(live_config):
    c,base=live_config
    assert '9 beams' in validate_data(c,base)
    flow=Workflow();flow.accept(0,c,base_directory=base)
    descriptor=flow.preparation_descriptor
    replay=apply_preparation(c,descriptor)
    assert replay==c
    for name in ('intensities','variances','weights'):
        np.testing.assert_array_equal(getattr(read_data(replay,base),name),getattr(read_data(c,base),name))
    original=descriptor.read_bytes()
    record=json.loads(original)
    metadata=json.loads(Path(record['preparation_metadata']['path']).read_text())
    assert metadata['excluded_indices'] and metadata['zone_axis_selection']['confirmed']
    manifest=json.loads(Path(metadata['overview']['manifest']).read_text())
    assert manifest['live_preparation']['noise_source']=='Poisson / intensity fallback'
    assert manifest['live_preparation']['lattice_detection']['seed']['basis_vector_2_px']==[-9.,82.]
    assert len(manifest['rois'])==9
    assert record['scan']['corner_a']==list(c.scan.corner_a)
    model=atomic_model(c)
    flow.sync(model);assert flow.receipts[0] is not None and flow.snapshot==c
    compiled=compile_config_from_sources(model)
    indexed=compiled.crystallography.indexed_beams
    assert [(b.experimental_index,b.original_hkl) for b in indexed]==[
        (b.lattice_coordinates,(2*b.lattice_coordinates[0],2*b.lattice_coordinates[1],0)) for b in c.crystallography.beams]
    assert descriptor.read_bytes()==original
    modified=replace(model,crystallography=replace(model.crystallography,beams=model.crystallography.beams[::-1]))
    flow.sync(modified);assert flow.snapshot is None and flow.receipts==[None]*3
    with pytest.raises(ValueError,match='identity/order'):validate_data(modified,base)


def test_live_seed_changes_a10_and_archived_corners_rejected(live_config):
    from qcbed.data_preparation.scan_conversion import convert_scan_geometry
    c,base=live_config;_,m=restored_overview(c,base)
    f=m['reciprocal_frame'];origin=np.asarray(f['zone_axis_origin_px'])
    assert f['basis_vector_1_px']==[81.,9.] and f['basis_vector_2_px']==[-9.,82.]
    def scan(g2):
        return convert_scan_geometry(reference_point_px=f['central_disk_center_px'],zone_axis_origin_px=tuple(origin),
            reference_point_1_px=tuple(origin+f['basis_vector_1_px']),reference_point_2_px=tuple(origin+g2),
            selection_width=76,selection_height=76,final_binning=2).scan_geometry
    current=scan(np.array([-9.,82.]));old=scan(np.array([-10.,81.]))
    assert current.corner_a==c.scan.corner_a and current.corner_b==c.scan.corner_b and current.corner_c==c.scan.corner_c
    assert old.corner_a!=current.corner_a
    with pytest.raises(ValueError,match='live A10'):validate_data(replace(c,scan=replace(c.scan,corner_a=old.corner_a)),base)


def test_real_abc_gui_initialization(live_config,monkeypatch):
    pytest.importorskip('PySide6')
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    from qcbed.grouped_run import load_grouped_checkpoint,grouped_history
    c,base=live_config;app=QApplication.instance() or QApplication([])
    w=ConfigWindow();w.document=EditorDocument(c,base/'project.yml');w.refresh()
    try:
        assert w.validate_stage(0),w.stage_status[0].text()
        w.preparation.install(str(w.workflow.preparation_descriptor))
        assert w.document.config==c and w.preparation.mode.currentIndex()==1
        w.document.commit_structure(atomic_model(c));w.refresh()
        assert w.validate_stage(1),w.stage_status[1].text()
        assert w.validate_c_part('simulation'),w.stage_status[2].text()
        assert w.validate_stage(2),w.stage_status[2].text()
        assert not w.run_action.isEnabled()
        assert w.validate_output_directory()
        w.check_run_readiness()
        assert w.run_action.isEnabled()
        assert not w.resume_action.isEnabled()
        from qcbed.gui.output_paths import modern_execution_config
        snapshot=prepare_run_snapshot(modern_execution_config(w.document.config,base),base_directory=base)
        run_directory=Path(snapshot.compiled.configuration.run.output_directory)
        for name in ('intensities','variances','weights'):
            np.testing.assert_array_equal(getattr(snapshot.experimental,name),getattr(read_data(c,base),name))
        class StopAtStage(Exception):pass
        def stop(*a,**k):raise StopAtStage()
        monkeypatch.setattr('qcbed.grouped_run._run_stage',stop)
        with pytest.raises(StopAtStage):run_snapshot(snapshot)
        cp=load_grouped_checkpoint(run_directory,snapshot)
        assert cp.get('evaluations',0)==0 and grouped_history(run_directory,cp)==[]
        assert not (run_directory/'.run.lock').exists()
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close();w.deleteLater();app.processEvents()


def test_changed_live_seed_flows_through_preparation(monkeypatch,live_config):
    import qcbed.data_preparation.live as live
    from qcbed.gui.preparation_model import load_image
    original=live._live_arc_seed
    def changed(*args):
        seed=original(*args)
        return replace(seed,basis_vector_2_px=(-8.,82.))
    monkeypatch.setattr(live,'_live_arc_seed',changed)
    r=live.prepare_from_raw_cbed(load_image(ROOT/'reference_data/pattern.dat'),
        LiveSelection(DiskGeometry(38,38),(111,105,155,151),76,76,final_binning=2,
                      signal_mode=SignalMode.GRADIENT,symmetry_crop_half_size_px=(19,19)))
    c,_=live_config
    assert r.dataset.reciprocal_frame.basis_vector_2_px==(-8.,82.)
    assert r.dataset.scan_geometry.corner_a!=c.scan.corner_a
    assert r.dataset.scan_geometry.corner_b!=c.scan.corner_b


def test_replay_rejects_changed_provenance(live_config,tmp_path):
    c,base=live_config;flow=Workflow();flow.accept(0,c,base_directory=base)
    descriptor=json.loads(flow.preparation_descriptor.read_text())
    descriptor['preparation_metadata']['sha256']='0'*64
    bad=tmp_path/'changed.prepared.json';bad.write_text(json.dumps(descriptor))
    with pytest.raises(ValueError,match='checksum mismatch'):apply_preparation(c,bad)

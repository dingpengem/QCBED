"""Headless tests of mouse-independent A8 selection/provenance/replay."""
import csv
import json
from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from qcbed.data_preparation.origin_selection import *
from qcbed.data_preparation.origin import measure_symmetry_origin


@pytest.fixture
def setup(tmp_path):
    image=np.arange(12000,dtype=float).reshape(100,120)
    rng=np.random.default_rng(7);crop=rng.random((20,20));crop+=crop[::-1,::-1]
    image[40:60,55:75]=crop
    pattern=tmp_path/'pattern.dat';np.savetxt(pattern,image)
    catalogue=tmp_path/'catalogue.csv'
    with catalogue.open('w') as f:
        w=csv.writer(f);w.writerow(['beam_id','n1','n2','predicted_x','predicted_y','status'])
        w.writerows([['central',0,0,30,50,'reference'],['weak',0,1,90,50,'weak'],['chosen',1,0,60,50,'observed']])
    return pattern,catalogue,image,load_selection_catalogue(catalogue)[2],DiskGeometry(20,20)


def test_nearest_valid_not_mouse_center(setup):
    _,catalogue,_,beam,_=setup
    result,distance=nearest_beam(load_selection_catalogue(catalogue),(62,51),selection_radius_px=5)
    assert result==beam and result.center_px==(60,50) and distance==pytest.approx(np.sqrt(5))
    with pytest.raises(ValueError,match='within'): nearest_beam(load_selection_catalogue(catalogue),(90,50),selection_radius_px=5)
    result,_=nearest_beam(load_selection_catalogue(catalogue),(90,50),selection_radius_px=5,allow_unreliable=True)
    assert result.status=='weak'
    direct,_=nearest_beam(load_selection_catalogue(catalogue),(30,50),selection_radius_px=5)
    assert direct.beam_id=='central'


def test_roi_round_trip_and_follow_beam(setup):
    _,_,_,beam,_=setup
    bounds=(60.,45.,70.,55.)
    local=detector_to_local_roi(bounds,beam.center_px)
    assert local==(0,-5,10,5)
    assert local_to_detector_roi(local,beam.center_px)==bounds
    assert local_to_detector_roi(local,(62,53))==(62,48,72,58)
    assert snap_roi_to_grid((60.2,45.3,69.8,55.1))==bounds


@pytest.mark.parametrize('bounds',[(60,45,60,55),(-1,2,3,4),(30,5,95,80),(60.2,45,70,55)])
def test_invalid_roi(setup,bounds):
    _,_,image,beam,geometry=setup
    with pytest.raises(ValueError): prepare_roi(bounds,beam,geometry,image.shape)


def test_crop_margin_and_advanced_containment(setup):
    _,_,image,beam,geometry=setup
    with pytest.raises(ValueError,match='within'):
        prepare_roi((60,45,70,55),beam,geometry,image.shape,search_half_size_px=(6,5))
    with pytest.raises(ValueError,match='room'):
        prepare_roi((5,45,15,55),replace(beam,center_px=(10,50)),geometry,image.shape)


def test_record_replay_exact_backend_and_hash_guard(setup,tmp_path,monkeypatch):
    pattern,catalogue,image,beam,geometry=setup
    before=image.copy()
    record=make_selection(pattern,catalogue,beam,(62,51),(60,45,70,55),geometry,(30,50),selection_radius_px=5,
                          raw_roi_detector=(60.2,45.3,69.8,55.1))
    reference,config=prepare_roi(record['roi_detector'],beam,geometry,image.shape)
    direct=measure_symmetry_origin(image,geometry,(30,50),reference,config=config)
    first=execute_selection(record,tmp_path/'first')
    assert first==direct
    saved=tmp_path/'first'/'a8_selection.json'
    data=json.loads(saved.read_text())
    assert data['roi_local']==[0,-5,10,5]
    assert data['result']['zone_axis_origin_px']==list(direct.zone_axis_origin_px)
    monkeypatch.chdir(tmp_path)
    assert replay_selection(saved,tmp_path/'replayed')==direct
    np.testing.assert_array_equal(image,before)
    with pytest.raises(FileExistsError): replay_selection(saved,tmp_path/'replayed')
    pattern.write_text(pattern.read_text()+'\n')
    with pytest.raises(ValueError,match='content changed'): replay_selection(saved,tmp_path/'changed')
    assert not (tmp_path/'changed').exists()


def test_provenance_geometry_tamper_rejected(setup,tmp_path):
    pattern,catalogue,_,beam,geometry=setup
    record=make_selection(pattern,catalogue,beam,(60,50),(60,45,70,55),geometry,(30,50),selection_radius_px=5)
    record['roi_local']=[1,2,3,4]
    with pytest.raises(ValueError,match='provenance'): execute_selection(record,tmp_path/'bad')
    assert not (tmp_path/'bad').exists()


def test_interactive_script_constructs_without_selecting_real_data(setup,tmp_path,monkeypatch):
    import importlib.util
    import matplotlib.pyplot as plt
    from types import SimpleNamespace
    path=Path(__file__).resolve().parents[1]/'examples/a8_origin_interactive.py'
    spec=importlib.util.spec_from_file_location('a8_interactive',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    pattern,catalogue,_,_,_=setup
    monkeypatch.setattr(plt,'show',lambda *a,**k:None)
    module.interactive(SimpleNamespace(pattern=pattern,catalogue=catalogue,radii=(20,20),central_center=(30,50),
        selection_radius=20,allow_unreliable=False,search_half_size=None,crop_half_size=None,output=tmp_path/'unused'))
    assert not (tmp_path/'unused').exists()
    plt.close('all')


def test_replay_cli_without_interactive_backend(setup,tmp_path):
    import subprocess,sys,os
    pattern,catalogue,_,beam,geometry=setup
    record=make_selection(pattern,catalogue,beam,(60,50),(60,45,70,55),geometry,(30,50),selection_radius_px=5)
    execute_selection(record,tmp_path/'original')
    script=Path(__file__).resolve().parents[1]/'examples/a8_origin_interactive.py'
    completed=subprocess.run([sys.executable,str(script),'--selection',str(tmp_path/'original/a8_selection.json'),
        '--output',str(tmp_path/'cli')],cwd=tmp_path,env=dict(os.environ,MPLBACKEND='Agg'),capture_output=True,text=True)
    assert completed.returncode==0,completed.stderr
    assert 'Zone-axis origin' in completed.stdout
    first=json.loads((tmp_path/'original/a8_selection.json').read_text())
    replay=json.loads((tmp_path/'cli/a8_selection.json').read_text())
    assert replay==first
    assert (tmp_path/'cli/selection_context.png').is_file()


def test_direct_reference_roi_replay_and_catalogue_order(setup,tmp_path):
    from qcbed.data_preparation.origin import select_symmetry_reference
    from types import SimpleNamespace as NS
    pattern,catalogue,image,_,geometry=setup
    beams=load_selection_catalogue(catalogue)
    for order in (beams,beams[::-1]):
        beam,_=nearest_beam(order,(31,50),selection_radius_px=5)
        assert beam.lattice_indices==(0,0)
    adapter=NS(sites=tuple(NS(beam_id=b.beam_id,integer_indices=b.lattice_indices,
        predicted_center_px=b.center_px,status=b.status) for b in beams[::-1]))
    ref=select_symmetry_reference(adapter,beam_id='central',feature_offset_px=(0,0))
    assert ref.beam_center_px==(30,50)
    record=make_selection(pattern,catalogue,beam,(31,50),(28,48,32,52),geometry,(30,50),selection_radius_px=5)
    result=execute_selection(record,tmp_path/'direct')
    assert result.reference_kind==record['reference_kind']=='direct'
    assert result.zone_axis_origin_px==result.symmetry_center_px
    assert replay_selection(tmp_path/'direct/a8_selection.json',tmp_path/'replay_direct')==result
    saved=json.loads((tmp_path/'direct/a8_selection.json').read_text())
    assert saved['reference_kind']=='direct' and saved['beam']['lattice_indices']==[0,0]
    assert saved['roi_local']==[-2,-2,2,2]


@pytest.mark.parametrize('confirm_snap_first',[False,True])
def test_mouse_roi_previews_snap_without_applying_until_confirmation(setup,tmp_path,monkeypatch,confirm_snap_first):
    import importlib.util
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button,RectangleSelector
    from matplotlib.backend_bases import FigureCanvasBase
    from types import SimpleNamespace as NS
    path=Path(__file__).resolve().parents[1]/'examples/a8_origin_interactive.py'
    spec=importlib.util.spec_from_file_location('a8_preview',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    callbacks={}
    original_connect=FigureCanvasBase.mpl_connect
    def connect(self,event,func):
        if getattr(func,'__name__','')=='click': callbacks['beam']=func
        return original_connect(self,event,func)
    monkeypatch.setattr(FigureCanvasBase,'mpl_connect',connect)
    original_button=Button.on_clicked
    def button(self,func):
        callbacks[self.label.get_text()]=func
        return original_button(self,func)
    monkeypatch.setattr(Button,'on_clicked',button)
    original_init=RectangleSelector.__init__
    def selector(self,ax,onselect,*args,**kwargs):
        callbacks['rectangle']=onselect;callbacks['axes']=ax
        original_init(self,ax,onselect,*args,**kwargs)
    monkeypatch.setattr(RectangleSelector,'__init__',selector)
    monkeypatch.setattr(plt,'show',lambda *a,**k:None)
    pattern,catalogue,_,_,_=setup
    output=tmp_path/'accepted'
    module.interactive(NS(pattern=pattern,catalogue=catalogue,radii=(20,20),central_center=(30,50),selection_radius=20,
        allow_unreliable=False,search_half_size=None,crop_half_size=None,output=output))
    callbacks['beam'](NS(inaxes=callbacks['axes'],button=1,xdata=60,ydata=50))
    callbacks['rectangle'](NS(xdata=60.1,ydata=45.1),NS(xdata=60.2,ydata=45.2))
    callbacks['Accept'](None)
    assert not output.exists()  # collapsed snap cannot be accepted as a raw ROI
    callbacks['Reselect ROI'](None)
    callbacks['rectangle'](NS(xdata=60.2,ydata=45.3),NS(xdata=69.8,ydata=55.1))
    fig=callbacks['axes'].figure
    assert any('Preview ROI=(60.0, 45.0, 70.0, 55.0)' in t.get_text() for t in fig.texts)
    assert not output.exists()  # preview is not numerical acceptance
    assert any(p.get_linestyle()=='--' for p in callbacks['axes'].patches)
    if confirm_snap_first:
        callbacks['Snap ROI to grid'](None)
        assert not output.exists()
    callbacks['Accept'](None)
    saved=json.loads((output/'a8_selection.json').read_text())
    assert saved['raw_roi_detector']==[60.2,45.3,69.8,55.1]
    assert saved['roi_detector']==[60,45,70,55]
    assert saved['roi_local']==[0,-5,10,5]
    assert saved['result']['symmetry_center_px']==[65,50]
    plt.close('all')

"""A8 mouse selection with explicit pixel-grid confirmation and headless JSON replay."""
import argparse
import hashlib
from pathlib import Path
from qcbed.data_preparation.models import DiskGeometry
from qcbed.data_preparation.diagnostics import load_pattern
from qcbed.data_preparation.origin_selection import (load_selection_catalogue,nearest_beam,
    detector_to_local_roi,snap_roi_to_grid,prepare_roi,make_selection,execute_selection,replay_selection)


def interactive(args):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import RectangleSelector,Button,RangeSlider
    from matplotlib.patches import Ellipse,Rectangle
    pattern=load_pattern(args.pattern);beams=load_selection_catalogue(args.catalogue)
    input_hashes=tuple(hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (args.pattern,args.catalogue))
    geometry=DiskGeometry(*args.radii)
    fig,ax=plt.subplots(figsize=(11,9));fig.subplots_adjust(bottom=.31,top=.92)
    low,high=np.percentile(pattern,[1,99.5])
    displayed=ax.imshow(pattern,origin='upper',cmap='gray',vmin=low,vmax=high)
    ax.set(xlabel='x (pixel)',ylabel='y (pixel)',title='A8: click a suitable direct or diffracted beam, then drag its local feature ROI')
    message=fig.text(.08,.18,'Click a suitable direct or diffracted beam. No beam or feature is chosen automatically.',fontsize=9,va='top')
    state=dict(phase='beam',beam=None,roi=None,raw=None,click=None,artists=[],proposal=None,invalid_proposal=False)
    def tell(text): message.set_text(text);fig.canvas.draw_idle()
    def draw():
        for item in state['artists']: item.remove()
        state['artists']=[]
        beam=state['beam']
        if beam is not None:
            ellipse=Ellipse(beam.center_px,2*geometry.radius_x_px,2*geometry.radius_y_px,fill=False,color='cyan');ax.add_patch(ellipse)
            state['artists'] += [ellipse,ax.scatter(*beam.center_px,c='cyan',marker='+')]
        if state['roi'] is not None:
            x0,y0,x1,y1=state['roi'];patch=Rectangle((x0,y0),x1-x0,y1-y0,fill=False,color='orange');ax.add_patch(patch)
            state['artists'] += [patch,ax.scatter((x0+x1)/2,(y0+y1)/2,c='orange',marker='x')]
        if state['proposal'] is not None and state['proposal']!=state['roi']:
            x0,y0,x1,y1=state['proposal']
            patch=Rectangle((x0,y0),x1-x0,y1-y0,fill=False,color='lime',linestyle='--')
            ax.add_patch(patch);state['artists'].append(patch)
        fig.canvas.draw_idle()
    def preview():
        draw();b=state['beam'];r=state['proposal'] or state['roi'];local=detector_to_local_roi(r,b.center_px)
        text=f'{b.beam_id} {b.lattice_indices}; center={b.center_px}; status={b.status}\nRaw ROI={state["raw"]} (orange)\nPreview ROI={r}; local={local} (green dashed when pending)\nSize=({r[2]-r[0]:g}, {r[3]-r[1]:g}); center=({(r[0]+r[2])/2:g}, {(r[1]+r[3])/2:g})'
        try:
            reference,config=prepare_roi(r,b,geometry,pattern.shape,search_half_size_px=args.search_half_size,crop_half_size_px=args.crop_half_size)
            text+=f'\nOffset={reference.feature_offset_px}; valid preview. Accept confirms these bounds and runs A8; Snap confirms without running.'
        except ValueError as exc: text+='\n'+str(exc)
        tell(text)
    def rectangle(start,end):
        if state['phase']!='roi' or None in (start.xdata,start.ydata,end.xdata,end.ydata): return
        r=(min(start.xdata,end.xdata),min(start.ydata,end.ydata),max(start.xdata,end.xdata),max(start.ydata,end.ydata))
        if r[0]==r[2] or r[1]==r[3]: tell('Zero-size ROI: drag a rectangle.');return
        state.update(roi=r,raw=r,proposal=None,invalid_proposal=False,phase='preview');selector.set_active(False)
        try: state['proposal']=snap_roi_to_grid(r)
        except ValueError:
            state['invalid_proposal']=True
            draw();tell(f'Raw ROI={r}\nNearest integer ROI has zero size; reselect ROI.');return
        preview()
    selector=RectangleSelector(ax,rectangle,useblit=True,button=[1],interactive=False)
    selector.set_active(False)
    def click(event):
        if state['phase']!='beam' or event.inaxes!=ax or event.button!=1 or event.xdata is None: return
        if fig.canvas.manager.toolbar is not None and fig.canvas.manager.toolbar.mode: return
        try:
            point=(event.xdata,event.ydata)
            beam,distance=nearest_beam(beams,point,selection_radius_px=args.selection_radius,allow_unreliable=args.allow_unreliable)
            state.update(beam=beam,click=point,phase='roi',roi=None,raw=None,proposal=None);draw()
            # Enable after release so the selecting click cannot also start an ROI.
            tell(f'{beam.beam_id} {beam.lattice_indices}; center={beam.center_px}; distance={distance:.3f} px; {beam.status}\nRelease, then drag a rectangle around the LOCAL symmetry feature.')
        except ValueError as exc: tell(str(exc))
    def release(event):
        if state['phase']=='roi': selector.set_active(True)
    def reselect_roi(event):
        if state['phase']=='done': return
        if state['beam'] is None: tell('Select a beam first.');return
        state.update(phase='roi',roi=None,raw=None,proposal=None);selector.clear();selector.set_active(True);draw();tell('Drag a new feature rectangle.')
    def reselect_beam(event):
        if state['phase']=='done': return
        state.update(phase='beam',beam=None,roi=None,raw=None,proposal=None);selector.clear();selector.set_active(False);draw();tell('Click a suitable direct or diffracted beam.')
    def snap(event):
        if state['phase']!='preview': return
        try:
            state['roi']=state['proposal'] or snap_roi_to_grid(state['roi']);state['proposal']=None;selector.clear();preview()
        except ValueError as exc: tell(str(exc))
    def accept(event):
        if state['phase']!='preview': tell('Select a beam and draw a valid ROI first.');return
        try:
            if state['invalid_proposal']: raise ValueError('Nearest integer ROI has zero size; reselect ROI')
            # Accept is explicit confirmation of the displayed proposed bounds.
            accepted_roi=state['proposal'] or state['roi']
            record=make_selection(args.pattern,args.catalogue,state['beam'],state['click'],accepted_roi,geometry,args.central_center,
                selection_radius_px=args.selection_radius,allow_unreliable=args.allow_unreliable,raw_roi_detector=state['raw'],
                search_half_size_px=args.search_half_size,crop_half_size_px=args.crop_half_size)
            if (record['pattern_sha256'],record['catalogue_sha256'])!=input_hashes:
                raise ValueError('Displayed pattern/catalogue changed on disk; restart selection to review current inputs')
            result=execute_selection(record,args.output)
        except (ValueError,TypeError,OSError,FloatingPointError) as exc: tell(f'Not accepted: {exc}\nReselect ROI/beam or choose a fresh output directory.');return
        state.update(phase='done',roi=accepted_roi,proposal=None);selector.clear();selector.set_active(False);draw()
        ax.scatter(*result.symmetry_center_px,c='lime',marker='+',label='A8 symmetry')
        ax.scatter(*result.zone_axis_origin_px,c='magenta',marker='*',label='Derived origin')
        ax.scatter(*args.central_center,c='white',marker='+',label='Central disk');ax.legend()
        tell(f'Saved {args.output}\n{result.reference.beam_id}: symmetry={result.symmetry_center_px}; mismatch={result.mismatch:g}\nOrigin={result.zone_axis_origin_px}; ROI local={record["roi_local"]}\nTies={result.tied_minimum_count}; boundary={result.boundary_limited}; contrast={result.crop_has_contrast}')
        fig.savefig(args.output/'interactive_preview.png',dpi=150)
        result_fig,result_ax=plt.subplots(figsize=(13,8));result_ax.imshow(plt.imread(args.output/'symmetry_diagnostic.png'));result_ax.axis('off');result_fig.tight_layout();plt.show(block=False)
        print(message.get_text())
    fig.canvas.mpl_connect('button_press_event',click);fig.canvas.mpl_connect('button_release_event',release)
    buttons=[]
    for i,(label,callback) in enumerate((('Accept',accept),('Snap ROI to grid',snap),('Reselect ROI',reselect_roi),('Reselect beam',reselect_beam),('Cancel / Close',lambda e:plt.close(fig)))):
        button=Button(fig.add_axes([.06+i*.185,.03,.175,.045]),label);button.on_clicked(callback);buttons.append(button)
    contrast=RangeSlider(fig.add_axes([.16,.24,.67,.025]),'Display %',0,100,valinit=(1,99.5))
    def adjust(values):
        if values[0]<values[1]: displayed.set_clim(*np.percentile(pattern,values));fig.canvas.draw_idle()
    contrast.on_changed(adjust)
    plt.show()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,help='Replay a saved a8_selection.json without mouse input')
    p.add_argument('--pattern',type=Path,default=Path(__file__).resolve().parents[1]/'reference_data/pattern.dat')
    p.add_argument('--catalogue',type=Path)
    p.add_argument('--central-center',type=float,nargs=2)
    p.add_argument('--radii',type=float,nargs=2)
    p.add_argument('--selection-radius',type=float,default=38.)
    p.add_argument('--allow-unreliable',action='store_true',help='Explicitly allow weak/missing/edge-clipped catalogue sites')
    p.add_argument('--search-half-size',type=float,nargs=2)
    p.add_argument('--crop-half-size',type=float,nargs=2)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.selection:
        import matplotlib
        matplotlib.use('Agg')
        result=replay_selection(a.selection,a.output);print('Zone-axis origin:',result.zone_axis_origin_px);return
    if a.catalogue is None or a.central_center is None or a.radii is None:
        p.error('Interactive mode requires --catalogue, --central-center and --radii')
    import matplotlib
    if matplotlib.get_backend().lower() in {'agg','pdf','svg','ps','pgf','cairo','template'}:
        p.error('An interactive matplotlib backend is required; use MPLBACKEND=QtAgg with the optional GUI environment, or a native interactive backend')
    interactive(a)


if __name__=='__main__': main()

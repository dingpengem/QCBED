"""Read-only plots of existing globally normalized objective arrays."""
from pathlib import Path
import numpy as np


def masked_comparison(arrays, rows, columns):
    mask=np.asarray(arrays['weights'])!=0
    if mask.ndim!=2 or mask.shape[0]!=rows*columns:
        raise ValueError('scan dimensions do not match tilt/beam arrays')
    result={'mask':mask.T.reshape(-1,rows,columns).copy()}
    for name in ('experimental','scaled_calculated','residual'):
        a=np.asarray(arrays[name])
        if a.shape!=mask.shape:raise ValueError('comparison array shape mismatch')
        result[name]=np.where(mask,a,np.nan).T.reshape(-1,rows,columns).copy()
    return result


def save_masked_pattern_match(arrays,rows,columns,labels,destination,*,title,ynorm):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data=masked_comparison(arrays,rows,columns)
    if len(labels)!=len(data['mask']):raise ValueError('beam labels do not match arrays')
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(destination.with_suffix('.npz'),**data,ynorm=ynorm)
    # Common display scale across all beams; no per-beam scientific scaling.
    intensity=np.concatenate([data[n][data['mask']] for n in ('experimental','scaled_calculated')])
    low,high=float(intensity.min()),float(intensity.max())
    limit=float(np.abs(data['residual'][data['mask']]).max()) or 1.
    fig,axes=plt.subplots(len(labels),3,figsize=(9,2*len(labels)),squeeze=False,layout='constrained')
    for b,label in enumerate(labels):
        for j,name in enumerate(('experimental','scaled_calculated','residual')):
            cmap=plt.get_cmap('coolwarm' if j==2 else 'viridis').with_extremes(bad='#dddddd')
            axes[b,j].imshow(np.ma.masked_invalid(data[name][b]),origin='upper',interpolation='nearest',
                cmap=cmap,vmin=-limit if j==2 else low,vmax=limit if j==2 else high)
            axes[b,j].set_xticks([]);axes[b,j].set_yticks([])
            if b==0:axes[b,j].set_title(('Experimental','Calculated × global YNORM','Calculated − experimental')[j])
        axes[b,0].set_ylabel(label)
    fig.suptitle(f'{title}; global YNORM={ynorm:.12g}\nGray = excluded support; shared intensity/residual scales')
    fig.savefig(destination,dpi=140);plt.close(fig)
    return data

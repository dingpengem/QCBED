"""A7 geometric enumeration and visible physical circumference (no image sampling)."""
from itertools import product
import numpy as np
from ..detector.models import _integer


def boundary_fraction(center, shape, geometry, *, margin=0.):
    """Arc-length-weighted midpoint quadrature, 4096 phases on the nominal ellipse.

    Pixel-center frame is [0,nx-1] × [0,ny-1]. A margin excludes derivative
    stencils near detector borders. No padded pixels contribute evidence.
    """
    t=(np.arange(4096)+.5)*(2*np.pi/4096)
    rx,ry=geometry.radius_x_px,geometry.radius_y_px
    x=center[0]+rx*np.cos(t);y=center[1]+ry*np.sin(t)
    speed=np.hypot(rx*np.sin(t),ry*np.cos(t))
    valid=(x>=margin)&(x<=shape[1]-1-margin)&(y>=margin)&(y<=shape[0]-1-margin)
    return float(speed[valid].sum()/speed.sum())


def enumerate_boundary_sites(origin, matrix, shape, geometry, *, max_sites=100000):
    """All nominal ellipses with nonzero in-frame circumference, even outside centers."""
    shape=tuple(_integer(n,'image dimension') for n in shape)
    if len(shape)!=2 or min(shape)<1: raise ValueError('Expected positive (ny,nx) shape')
    matrix=np.asarray(matrix,dtype=float);origin=np.asarray(origin,dtype=float)
    if matrix.shape!=(2,2) or not np.isfinite(matrix).all() or np.linalg.cond(matrix)>1e6:
        raise ValueError('Invalid lattice basis')
    radii=np.array([geometry.radius_x_px,geometry.radius_y_px])
    lower=-radii;upper=np.array(shape[::-1])-1+radii
    corners=np.array(list(product((lower[0],upper[0]),(lower[1],upper[1]))))
    indices=np.linalg.solve(matrix,(corners-origin).T).T
    lo=np.floor(indices.min(axis=0)).astype(int);hi=np.ceil(indices.max(axis=0)).astype(int)
    if int(hi[0]-lo[0]+1)*int(hi[1]-lo[1]+1)>max_sites:
        raise ValueError('Enumeration exceeds max_enumeration_sites; no truncation')
    sites=[]
    for i in range(lo[0],hi[0]+1):
        for j in range(lo[1],hi[1]+1):
            point=origin+matrix@[i,j]
            if np.any(point<lower) or np.any(point>upper): continue
            if boundary_fraction(point,shape,geometry)>0:
                sites.append(((i,j),tuple(float(v) for v in point)))
    sites.sort(key=lambda s:(s[0]!=(0,0),float(np.sum((np.array(s[1])-origin)**2)),s[0]))
    return tuple(sites)

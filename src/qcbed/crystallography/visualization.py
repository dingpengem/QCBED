"""Optional Matplotlib 3-D diagnostic; never imported by lattice/structure physics."""
from itertools import product
import numpy as np
from .structure import AtomicStructureModel
from .replication import replica_offsets, replicated_cartesian_positions


def plot_atomic_structure(structure: AtomicStructureModel, repeats=(1,1,1), *, ax=None):
    """Plot cell edges and unwrapped Cartesian atoms (Å); return (figure, axes).

    Replicas translate by integer direct lattice vectors, without modifying sites.
    Draws every replicated unit-cell frame. Does not display or save implicitly.
    """
    if not isinstance(structure, AtomicStructureModel):
        raise TypeError('structure must be AtomicStructureModel')
    offsets = replica_offsets(repeats)
    import matplotlib.pyplot as plt
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        if getattr(ax, 'name', None) != '3d':
            raise ValueError('ax must be a Matplotlib 3-D axes')
        fig = ax.figure
    a = structure.direct_lattice.vectors
    positions = replicated_cartesian_positions(structure, offsets)
    repeated_symbols = [site.symbol for _ in offsets for site in structure.sites]
    corners = np.array(list(product((0,1), repeat=3)))
    extent = []
    for offset in offsets:
        cart = (corners+offset) @ a
        extent.extend(cart)
        for i,start in enumerate(corners):
            for axis in range(3):
                if start[axis] == 0:
                    end = start.copy(); end[axis] = 1
                    edge = np.array([start+offset,end+offset]) @ a
                    ax.plot(*edge.T, color='gray', linewidth=.7)
    symbols = tuple(dict.fromkeys(site.symbol for site in structure.sites))
    for i,symbol in enumerate(symbols):
        cart = positions[np.array([v == symbol for v in repeated_symbols])]
        extent.extend(cart)
        ax.scatter(*cart.T, label=symbol, color=plt.get_cmap('tab20')(i % 20),
                   marker=('o','^','s','D')[i//20 % 4], s=45)
    extent = np.array(extent)
    spans = np.ptp(extent, axis=0)
    ax.set_box_aspect(spans)
    for setter,axis in zip((ax.set_xlim,ax.set_ylim,ax.set_zlim),range(3)):
        setter(extent[:,axis].min()-.05*spans[axis],extent[:,axis].max()+.05*spans[axis])
    ax.set(xlabel='x (Å)', ylabel='y (Å)', zlabel='z (Å)')
    ax.legend(title='Element')
    return fig, ax


def plot_projected_atomic_structure(structure, orientation, repeats=(1,1,1), *, ax=None):
    """Unmerged 2-D beam coordinates in Å, equal aspect; return (figure,axes)."""
    from .projection import project_atomic_structure
    sites = project_atomic_structure(structure, orientation, repeats)
    import matplotlib.pyplot as plt
    if ax is None:
        fig,ax = plt.subplots()
    else:
        if getattr(ax, 'name', None) != 'rectilinear':
            raise ValueError('ax must be a 2-D rectilinear Matplotlib axes')
        fig = ax.figure
    for i,symbol in enumerate(dict.fromkeys(s.symbol for s in sites)):
        rows = [s for s in sites if s.symbol == symbol]
        ax.scatter([s.x for s in rows],[s.y for s in rows],label=symbol,
            color=plt.get_cmap('tab20')(i % 20),marker=('o','^','s','D')[i//20 % 4])
    zone = ' '.join(str(v) for v in orientation.zone_axis.uvw)
    ax.set(xlabel='x_beam (Å)',ylabel='y_beam (Å)',title=f'Zone axis [{zone}]')
    ax.set_aspect('equal',adjustable='box');ax.legend(title='Element')
    return fig,ax

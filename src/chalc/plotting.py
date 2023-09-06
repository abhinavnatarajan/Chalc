from __future__ import annotations

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import seaborn as sns
from numbers import Real
from .sixpack import DiagramEnsemble, _num_colours_in_bitmask
from .filtration import FilteredComplex

plt.rcParams["animation.html"] = "jshtml"

# __all__ = ['plot_sixpack', 'animate_filtration']

__doc__ = 'Plotting and visualisation utilities.'

def plot_sixpack(dgms : DiagramEnsemble,
                 truncation : Real, 
                 max_diagram_dim : int = 2) -> None :
    """
    Plots the 6-pack of persistence diagrams returned by :func:`compute <.sixpack.compute>`.

    Parameters
    ----------
    dgms : 
        The persistence diagrams.
    truncation :
        The maximum entrance time for which the diagrams are plotted.
    max_diagram_dim :
        Maximum homological dimension for which the diagrams are plotted.

    """

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[3.5 * 3, 3.5 * 2])
    
    def applied_plot(
        dgm, ax, title, dim_shift=0, max_dim=max_diagram_dim, legend=False
    ):
        _plot_diagram(
            dgm,
            dgms.entrance_times,
            dgms.dimensions,
            truncation,
            ax=ax,
            title=title,
            dim_shift=dim_shift,
            max_dim=max_dim,
            legend=legend,
        )

    applied_plot(
        dgms.ker,
        axes[0][0],
        "Kernel",
        dim_shift=1,
    )
    applied_plot(
        dgms.rel,
        axes[0][1],
        "Relative",
        max_dim=max_diagram_dim + 1,
    )
    applied_plot(
        dgms.cok,
        axes[0][2],
        "Cokernel",
    )
    applied_plot(
        dgms.cod,
        axes[1][0],
        "Domain",
    )
    applied_plot(
        dgms.im,
        axes[1][1],
        "Image",
    )
    applied_plot(
        dgms.dom,
        axes[1][2],
        "Codomain",
        legend=True,
    )
    return fig, axes


def _plot_diagram(
    diagram,
    entrance_times,
    dimensions,
    truncation,
    max_dim=2,
    ax=None,
    title=None,
    legend=True,
    dim_shift=0,
):
    # Truncate all times
    entrance_times = [
        et if et <= truncation * 1.05 else truncation * 1.05 for et in entrance_times
    ]
    all_pts = np.array(
        [
            [
                entrance_times[birth_idx],
                entrance_times[death_idx],
                dimensions[birth_idx] - dim_shift,
            ]
            for birth_idx, death_idx in diagram.paired
            if entrance_times[death_idx] != entrance_times[birth_idx]
            and dimensions[birth_idx] - dim_shift <= max_dim
        ]
        + [
            [
                entrance_times[birth_idx],
                truncation * 1.05,
                dimensions[birth_idx] - dim_shift,
            ]
            for birth_idx in diagram.unpaired
            if dimensions[birth_idx] - dim_shift <= max_dim
        ]
    )
    df = pd.DataFrame(data=list(all_pts), columns=["Birth", "Death", "Dimension"])
    ret_ax = sns.scatterplot(
        data=df, x="Birth", y="Death", hue="Dimension", ax=ax, legend=legend
    )
    ret_ax.set(xlabel=None)
    ret_ax.set(ylabel=None)
    ret_ax.set(title=title)
    if legend:
        sns.move_legend(ret_ax, "lower right")
    handle = ax if ax is not None else plt
    handle.plot([0, truncation * 1.05], [0, truncation * 1.05], "m", alpha=0.4)
    handle.plot(
        [0, 0, truncation * 1.05],
        [0, truncation * 1.05, truncation * 1.05],
        "m--",
        alpha=0.4,
    )

def animate_filtration(
        K : FilteredComplex, 
        points : numpy.ndarray[numpy.float64[m, n]],
        filtration_times : list[float],
        length : float):
    """
    Create animation of 2-skeleton of filtered simplicial complex.

    Parameters
    ----------
    K : 
        A filtered complex.
    points : 
        The vertices of ``K`` as a (2, N) numpy matrix.
    filtration_times :
        List of filtration times for which to draw animation frames.
    length :  
        Total length of the animation in seconds. 

    Returns
    -------
    matplotlib.animation.FuncAnimation
    """
    if len(points.shape) != 2:
        raise NotImplementedError
    
    fig, ax = plt.subplots()
    plot_colours = np.array(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    vertex_colours = []
    for idx, simplex in K.simplices[0].items():
        i = int(np.round(np.log2(simplex.colours)))
        vertex_colours.append(i)
    
    ax.scatter(points[0, :], points[1, :], c=plot_colours[vertex_colours], s=10)

    lines = dict()
    patches = dict()
    for idx, simplex in K.simplices[1].items():
        if _num_colours_in_bitmask(simplex.colours) == 1:
            i = int(np.round(np.log2(simplex.colours)))
            colour = plot_colours[i]
        else:
            colour = 'grey'
        lines[idx] = ax.plot([], [], c=colour, alpha=0.5, linewidth=1)[0]

    for idx, simplex in K.simplices[2].items():
        if _num_colours_in_bitmask(simplex.colours) == 1:
            i = int(np.round(np.log2(simplex.colours)))
            colour = plot_colours[i]
        else:
            colour = 'grey'
        patches[idx] = ax.fill([], [], c=colour, alpha=0.2)[0]

    def update(time):
        for idx, simplex in K.simplices[1].items():
            if simplex.filtration_value <= time:
                lines[idx].set_xdata(points[0, simplex.vertices])
                lines[idx].set_ydata(points[1, simplex.vertices])

        for idx, simplex in K.simplices[2].items():
            if simplex.filtration_value <= time:
                patches[idx].set_xy(points[:, simplex.vertices].T)

    interval = int(np.round(len(filtration_times) / length))
    return animation.FuncAnimation(
        fig=fig, func=update, 
        frames=filtration_times, 
        interval=interval)
    
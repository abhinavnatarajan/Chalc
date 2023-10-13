from __future__ import annotations

from matplotlib.figure import FigureBase
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from numpy.typing import NDArray
from typing import Annotated
import pandas as pd
import seaborn as sns
from .sixpack import DiagramEnsemble, _num_colours_in_bitmask, _bitmask_to_colours, _colours_are_subset, _colours_to_bitmask
from .filtration import FilteredComplex
from ._utils import interpolate_docstring
import scipy

plt.rcParams["animation.html"] = "jshtml"

__doc__ = 'Plotting and visualisation utilities.'

def plot_sixpack(dgms : DiagramEnsemble,
                 *,
                 truncation : float | None = None,
                 max_diagram_dimension : int | None = None) -> tuple[FigureBase, Annotated[NDArray, "Axes"]] :
    """
    Plots the 6-pack of persistence diagrams returned by :func:`compute <.sixpack.compute>`.

    Args:
        dgms : The six-pack of persistence diagrams.

    Keyword Args:
        truncation : The maximum entrance time for which the diagrams are plotted. A sensible default will be calculated if not provided.
        max_diagram_dimension : The maximum homological dimension for which to plot points. If not provided, all dimensions will be included in the plots.
    """
    if truncation is None:
        truncation = _get_truncation(dgms.entrance_times)
    if max_diagram_dimension is None:
        max_diagram_dimension = max(dgms.dimensions)
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[3.5 * 3, 3.5 * 2])
    def applied_plot(
        dgm, ax, title, dim_shift=0, max_dim=max_diagram_dimension, legend=False, **kwargs
    ):
        _plot_diagram(
            dgm,
            dgms.entrance_times,
            dgms.dimensions,
            truncation,
            ax        = ax,
            title     = title,
            dim_shift = dim_shift,
            max_dim   = max_dim,
            legend    = legend,
            **kwargs
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
        max_dim = max_diagram_dimension + 1,
        legend  = True
    )
    applied_plot(
        dgms.cok,
        axes[0][2],
        "Cokernel",
    )
    applied_plot(
        dgms.dom,
        axes[1][0],
        "Domain",
    )
    applied_plot(
        dgms.im,
        axes[1][1],
        "Image",
    )
    applied_plot(
        dgms.cod,
        axes[1][2],
        "Codomain",
        legend = True,
    )
    return fig, axes

@interpolate_docstring
def plot_diagram(dgms: DiagramEnsemble,
                 diagram_name: str,
                 *,
                 truncation : float | None = None,
                 max_diagram_dimension : int | None = None ,
                 ) -> tuple[FigureBase, Annotated[NDArray, "Axes"]] :
    """
    Plot a specific diagram from a sixpack.

    Args:
        dgms : The six-pack of persistence diagrams.
        diagram_name : One of ``${str(DiagramEnsemble.diagram_names)}``.

    Keyword Args:
        truncation : The maximum entrance time for which the diagrams are plotted. A sensible default will be calculated if not provided.
        max_diagram_dimension : The maximum homological dimension for which to plot points. If not provided, all dimensions will be included in the plots.
    """
    if not diagram_name in DiagramEnsemble.diagram_names:
        raise KeyError("Invalid diagram name!")
    titles = {
        'ker': 'Kernel',
        'cok': 'Cokernel',
        'im' : 'Image',
        'dom': 'Domain',
        'cod': 'Codomain',
        'rel': 'Relative',
    }
    if truncation is None:
        truncation = _get_truncation(dgms.entrance_times)
    if max_diagram_dimension is None:
        max_diagram_dimension = max(dgms.dimensions)
    fig, ax = plt.subplots()
    _plot_diagram(
        getattr(dgms, diagram_name),
        dgms.entrance_times,
        dgms.dimensions,
        truncation,
        ax        = ax,
        dim_shift = 1 if diagram_name == 'ker' else 0,
        legend    = True,
        title     = titles[diagram_name],
        max_dim   = max_diagram_dimension
    )
    return fig, ax

def _plot_diagram(
    diagram,
    entrance_times,
    dimensions,
    truncation,
    max_dim   = 2,
    ax        = None,
    title     = None,
    legend    = True,
    dim_shift = 0,
    **kwargs
):
    plot_colours = plt.rcParams['axes.prop_cycle'].by_key()['color'][2 : ]
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
        data=df, x="Birth", y="Death", hue="Dimension", ax=ax, legend=legend, palette=plot_colours, **kwargs
    )
    ret_ax.set(xlabel=None)
    ret_ax.set(ylabel=None)
    ret_ax.set(title=title)
    ret_ax.set_aspect('equal')
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

def draw_filtration(
        K : FilteredComplex,
        points : NDArray[np.float64],
        time :  float,
        *,
        include_colours : list[int] | None = None) -> tuple[FigureBase, Axes] :
    """
    Visualise a filtration at given time, optionally including only certain colours.

    Parameters
    ----------
    K :
        A filtered complex.
    points :
        The vertices of ``K`` as a :math:`2\\times N` numpy matrix.
    time :
        Filtration times for which to draw simplices.

    Keyword Args
    ------------
    include_colours :
        Optional list of colours to include. If not specified then all colours will be drawn.
    """
    if len(points.shape) != 2:
        raise NotImplementedError

    if include_colours is None:
        include_colours = list(set([_bitmask_to_colours(vertex.colours)[0] for vertex in K.simplices[0].values()]))

    include_colours_bitmask = _colours_to_bitmask(include_colours)

    fig : FigureBase
    ax : Axes
    fig, ax = plt.subplots()
    plot_colours = np.array(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    # Plot the vertices
    vertex_colours = []
    vertices_to_plot = []
    for idx, vertex in K.simplices[0].items():
        if (_colours_are_subset(vertex.colours, include_colours_bitmask)) and vertex.filtration_value <= time:
            vertex_colours += _bitmask_to_colours(vertex.colours)
            vertices_to_plot.append(idx)

    vertex_colours = np.array(vertex_colours)
    ax.scatter(points[0, vertices_to_plot], points[1, vertices_to_plot], c=list(plot_colours[vertex_colours]), s=10)

    # Plot the edges
    for idx, simplex in K.simplices[1].items():
        if _colours_are_subset(simplex.colours, include_colours_bitmask) and simplex.filtration_value <= time:
            if _num_colours_in_bitmask(simplex.colours) == 1:
                colour = plot_colours[_bitmask_to_colours(simplex.colours)[0]]
                alpha = 0.5
            else:
                colour = 'black'
                alpha = 0.2
            ax.plot(points[0, simplex.vertices], points[1, simplex.vertices], c=colour, alpha=alpha, linewidth=1)

    for idx, simplex in K.simplices[2].items():
        if _colours_are_subset(simplex.colours, include_colours_bitmask) and simplex.filtration_value <= time:
            if _num_colours_in_bitmask(simplex.colours) == 1:
                colour = plot_colours[_bitmask_to_colours(simplex.colours)[0]]
            else:
                colour = 'grey'
            ax.fill(points[0, simplex.vertices], points[1, simplex.vertices], c=colour, alpha=0.2)

    ax.set_aspect('equal')
    # ax.set_xlabel('Time = ' + f"{0.0:.4f}")

    return fig, ax

def animate_filtration(
        K : FilteredComplex,
        points : NDArray[np.float64],
        *,
        filtration_times : list[float],
        animation_length : float) -> animation.FuncAnimation :
    """
    Create animation of 2-skeleton of filtered simplicial complex.

    Parameters
    ----------
    K :
        A filtered complex.
    points :
        The vertices of ``K`` as a :math:`2\\times N` numpy matrix.

    Keyword Args
    ------------
    filtration_times :
        List of filtration times for which to draw animation frames.
    animation_length :
        Total length of the animation in seconds.
    """
    if len(points.shape) != 2:
        raise NotImplementedError

    fig : FigureBase
    ax : Axes
    fig, ax = plt.subplots()
    plot_colours = np.array(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    vertex_colours = []
    for idx, simplex in K.simplices[0].items():
        i = _bitmask_to_colours(simplex.colours)[0]
        vertex_colours.append(i)

    ax.scatter(points[0, :], points[1, :], c=list(plot_colours[vertex_colours]), s=10)

    lines = dict()
    patches = dict()
    for idx, simplex in K.simplices[1].items():
        if _num_colours_in_bitmask(simplex.colours) == 1:
            i = _bitmask_to_colours(simplex.colours)[0]
            colour = plot_colours[i]
            alpha = 0.5
        else:
            colour = 'black'
            alpha = 0.2
        lines[idx] = ax.plot([], [], c=colour, alpha=alpha, linewidth=1)[0]

    for idx, simplex in K.simplices[2].items():
        if _num_colours_in_bitmask(simplex.colours) == 1:
            i = _bitmask_to_colours(simplex.colours)[0]
            colour = plot_colours[i]
        else:
            colour = 'grey'
        patches[idx] = ax.fill([], [], c=colour, alpha=0.2)[0]

    ax.set_aspect('equal')
    ax.set_xlabel('Time = ' + f"{0.0:.4f}")

    def update(time : float) -> None:
        for idx, simplex in K.simplices[1].items():
            if simplex.filtration_value <= time:
                lines[idx].set_xdata(points[0, simplex.vertices])
                lines[idx].set_ydata(points[1, simplex.vertices])

        for idx, simplex in K.simplices[2].items():
            if simplex.filtration_value <= time:
                patches[idx].set_xy(points[:, simplex.vertices].T)

        ax.set_xlabel(f'Time = {time:.4f}')

    interval = int(np.round(animation_length * 1000 / len(filtration_times)))
    return animation.FuncAnimation(
        fig      = fig,
        func     = update,
        frames   = filtration_times,
        interval = interval)

def _get_truncation(entrance_times : list[float]) :
    ub = scipy.stats.scoreatpercentile(entrance_times, 95)
    return ub if (m := max(entrance_times) > 1.5 * ub) else m

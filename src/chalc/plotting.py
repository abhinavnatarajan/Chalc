"""Plotting and visualisation utilities."""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence, Set
from itertools import product
from typing import Literal

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from pandas import DataFrame

from chalc.filtration import FilteredComplex

from .sixpack import (
	DiagramEnsemble,
	SimplexPairings,
	_bitmask_to_colours,
	_colours_are_subset,
	_colours_to_bitmask,
	_num_colours_in_bitmask,
)

plt.rcParams["animation.html"] = "jshtml"


def plot_sixpack(
	dgms: DiagramEnsemble,
	truncations: Mapping[DiagramEnsemble.DiagramName, float] | float | None = None,
	dimensions: Set[int] | int | None = None,
	tolerance: float = 0,
) -> tuple[Figure, np.ndarray[tuple[int, ...], np.dtype[Axes]]]:
	"""
	Plots the 6-pack of persistence diagrams returned by :func:`compute <.sixpack.compute>`.

	Args:
		dgms        : The 6-pack of persistence diagrams.
		truncations : The maximum entrance time upto which features are plotted. Can be either a single value, or a mapping from diagram names to values. A sensible default will be calculated if not provided.
		dimensions  : The homological dimensions for which to plot features. If not provided, all dimensions will be included in the plots.
		tolerance   : Only features with persistence greater than this value will be plotted.
	"""
	if dimensions is None:
		if dgms:
			dimensions = set(range(max(dgms.dimensions)))
			# The relative diagram may contain features in one dimension higher than that of the filtration.
			dims = {
				name: (dimensions if name != "rel" else dimensions | {max(dimensions) + 1})
				for name in DiagramEnsemble.diagram_names
			}
		else:
			dimensions = set()
	elif isinstance(dimensions, int):
		dimensions = {
			dimensions,
		}

	# In general, the dimension of a feature represented by a simplex pair (a, b)
	# is dim(a), unless the diagram is in the kernel, in which case
	# it is dim(a) - 1.
	dim_shift = {name: (1 if name == "ker" else 0) for name in DiagramEnsemble.diagram_names}

	if truncations is None:
		truncations = dict()
	if isinstance(truncations, Mapping):
		truncs = dict()
		for diagram_name in DiagramEnsemble.diagram_names:
			if diagram_name in truncations.keys() and isinstance(truncations[diagram_name], float):
				truncs[diagram_name] = truncations[diagram_name]
				continue
			ycoords = [
				dgms._entrance_times[death_simplex]
				for birth_simplex, death_simplex in dgms[diagram_name].paired
				if dgms._dimensions[birth_simplex] - dim_shift[diagram_name] in dims[diagram_name]
			]
			truncs[diagram_name] = _get_truncation(ycoords) if ycoords else 1.0
	elif isinstance(truncations, float):
		truncs = {diagram_name: truncations for diagram_name in DiagramEnsemble.diagram_names}

	plot_pos = {
		name: (val[1], val[0])
		for name, val in zip(DiagramEnsemble.diagram_names, product((0, 1, 2), (0, 1)))
	}
	plot_titles = {
		"ker": "Kernel",
		"cok": "Cokernel",
		"im": "Image",
		"dom": "Domain",
		"cod": "Codomain",
		"rel": "Relative",
	}
	points_legend = {
		name: (True if name == "rel" else False) for name in DiagramEnsemble.diagram_names
	}
	fig, axes = plt.subplots(
		nrows=2, ncols=3, figsize=[3.5 * 3, 3.5 * 2], sharex=False, sharey=False
	)
	for diagram_name in DiagramEnsemble.diagram_names:
		_plot_diagram(
			diagram=dgms[diagram_name],
			entrance_times=list(dgms._entrance_times),
			dimensions=list(dgms._dimensions),
			truncation=truncs[diagram_name],
			dims=dims[diagram_name],
			ax=axes[plot_pos[diagram_name]],
			dim_shift=dim_shift[diagram_name],
			title=plot_titles[diagram_name],
			points_legend=points_legend[diagram_name],
			lines_legend=False,
			tolerance=tolerance,
		)
	legends = dict()
	for ax in fig.axes:
		handles, labels = ax.get_legend_handles_labels()
		for handle, label in zip(handles, labels):
			legends[label] = handle
		ax.tick_params(labelleft=True, labelbottom=True)
	axes[plot_pos["rel"]].get_legend().remove()
	fig.legend(handles=legends.values(), loc="lower center", ncol=5)
	return fig, axes


def plot_diagram(
	dgms: DiagramEnsemble,
	diagram_name: DiagramEnsemble.DiagramName,
	truncation: float | None = None,
	dimensions: Set[int] | int | None = None,
	ax: Axes | None = None,
	tolerance: float = 0,
) -> Axes | None:
	"""
	Plot a specific diagram from a 6-pack.

	Args:
		dgms         : The 6-pack of persistence diagrams.
		diagram_name : One of ``'ker'``, ``'cok'``, ``'dom'``, ``'cod'``, ``'im'``, or ``'rel'``.
		truncation   : The maximum entrance time for which the diagrams are plotted. A sensible default will be calculated if not provided.
		dimensions   : The homological dimensions for which to plot features. If not provided, all dimensions will be included in the plots.
		ax           : A matplotlib axes object. If provided then the diagram will be plotted on the given axes.
		tolerance    : Only features with persistence greater than this value will be plotted.
	"""
	dgm = dgms[diagram_name]

	# In general, the dimension of a feature represented by a simplex pair (a, b)
	# is dim(a), unless the diagram is in the kernel, in which case
	# it is dim(a) - 1.
	dim_shift = 1 if diagram_name == "ker" else 0

	if dimensions is None:
		if dgm:
			dimensions = set(range(max(dgms._dimensions)))
		else:
			dimensions = set()
	elif isinstance(dimensions, int):
		dimensions = {
			dimensions,
		}

	if truncation is None:
		ycoord = [
			dgms._entrance_times[death_simplex]
			for birth_simplex, death_simplex in dgms[diagram_name].paired
			if dgms._dimensions[birth_simplex] - dim_shift in dimensions
		]
		truncation = _get_truncation(ycoord) if ycoord else 1.0

	if ax is None:
		ax1: Axes
		_, ax1 = plt.subplots()
	else:
		ax1 = ax

	titles: dict[DiagramEnsemble.DiagramName, str] = {
		"ker": "Kernel",
		"cok": "Cokernel",
		"im": "Image",
		"dom": "Domain",
		"cod": "Codomain",
		"rel": "Relative",
	}

	_plot_diagram(
		diagram=dgm,
		entrance_times=list(dgms._entrance_times),
		dimensions=list(dgms._dimensions),
		truncation=truncation,
		dims=dimensions,
		ax=ax1,
		title=titles[diagram_name],
		points_legend=True,
		lines_legend=True,
		dim_shift=dim_shift,
		tolerance=tolerance,
	)
	return ax1


def _plot_diagram(
	diagram: SimplexPairings,
	entrance_times: Sequence[float],  # entrance times of the simplices by index
	dimensions: Sequence[int],  # dimensions of the simplices by index
	truncation: float,  # max birth/death time to plot
	dims: Set[int],  # dimensions for which features are plotted
	ax: Axes,  # axes to plot on
	title: str,  # title of the plot
	points_legend: bool,  # show legend for points (which dimension)
	lines_legend: bool,  # show legend for lines (truncation, infinity)
	dim_shift: int,  # dimension shift if any; only relevant for kernel
	tolerance: float,
	**kwargs,
) -> None:
	# Truncate all times
	truncate_level = truncation * 1.04
	inf_level = truncation * 1.08
	entrance_times = [
		et if et <= truncation else truncate_level if et is not np.inf else inf_level
		for et in entrance_times
	]
	all_pts = [
		(
			entrance_times[birth_idx],
			entrance_times[death_idx],
			f"$H_{{{dimensions[birth_idx] - dim_shift}}}$",
		)
		for birth_idx, death_idx in diagram._paired
		if dimensions[birth_idx] - dim_shift in dims
		and entrance_times[death_idx] - entrance_times[birth_idx] > tolerance
	] + [
		(
			entrance_times[birth_idx],
			inf_level,
			f"$H_{{{dimensions[birth_idx] - dim_shift}}}$",
		)
		for birth_idx in diagram._unpaired
		if dimensions[birth_idx] - dim_shift in dims
	]
	df = DataFrame(data=all_pts, columns=["Birth", "Death", "Dimension"])
	df["Dimension"] = df["Dimension"].astype("category")
	ret_ax = sns.scatterplot(
		data=df, x="Birth", y="Death", hue="Dimension", ax=ax, legend=points_legend, **kwargs
	)
	ret_ax.set(xlabel=None)
	ret_ax.set(ylabel=None)
	ret_ax.set(title=title)
	ret_ax.set_aspect("equal")
	has_legend = any(True for c in ret_ax.get_children() if isinstance(c, Legend))
	if points_legend and has_legend:
		sns.move_legend(ret_ax, "lower right")
	handle = ax if ax is not None else plt
	handle.plot([0, inf_level], [0, inf_level], "k-", alpha=0.4)
	handle.plot(
		[0, truncate_level], [truncate_level, truncate_level], "k:", alpha=0.4, label="Truncated"
	)
	handle.plot([0, inf_level], [inf_level, inf_level], "k--", alpha=0.4, label="$\\infty$")
	handle.plot([0, 0], [0, inf_level], "k-", alpha=0.4)
	if lines_legend:
		handle.legend()


def draw_filtration(
	K: FilteredComplex,
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.float64]],
	time: float,
	include_colours: Collection[int] | None = None,
	ax: Axes | None = None,
) -> Axes:
	"""
	Visualise a 2D filtration at given time, optionally including only certain colours.

	Args:
		K               : A filtered 2-dimensional simplicial complex.
		points          : The vertices of ``K`` as a :math:`2\\times N` numpy matrix.
		time            : Filtration times for which to draw simplices.
		include_colours : Optional collection of colours to include. If not specified then all colours will be drawn.
		ax              : A matplotlib axes object. If provided then the diagram will be plotted on the given axes.
	"""
	if len(points.shape) != 2:
		raise NotImplementedError

	if include_colours is None:
		include_colours = set(
			[_bitmask_to_colours(vertex.colours)[0] for vertex in K.simplices[0].values()]
		)

	include_colours_bitmask = _colours_to_bitmask(include_colours)

	if ax is None:
		ax1: Axes
		_, ax1 = plt.subplots()
	else:
		ax1 = ax

	plot_colours = np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])

	# Plot the vertices
	vertex_colours = []
	vertices_to_plot = []
	for idx, vertex in K.simplices[0].items():
		if (
			_colours_are_subset(vertex.colours, include_colours_bitmask)
		) and vertex.filtration_value <= time:
			vertex_colours += _bitmask_to_colours(vertex.colours)
			vertices_to_plot.append(idx)

	vertex_colours = np.array(vertex_colours)
	ax1.scatter(
		points[0, vertices_to_plot],
		points[1, vertices_to_plot],
		c=list(plot_colours[vertex_colours]),
		s=10,
	)

	# Plot the edges
	for idx, simplex in K.simplices[1].items():
		if (
			_colours_are_subset(simplex.colours, include_colours_bitmask)
			and simplex.filtration_value <= time
		):
			if _num_colours_in_bitmask(simplex.colours) == 1:
				colour = plot_colours[_bitmask_to_colours(simplex.colours)[0]]
				alpha = 0.5
			else:
				colour = "black"
				alpha = 0.2
			ax1.plot(
				points[0, simplex.vertices],
				points[1, simplex.vertices],
				c=colour,
				alpha=alpha,
				linewidth=1,
			)

	for idx, simplex in K.simplices[2].items():
		if (
			_colours_are_subset(simplex.colours, include_colours_bitmask)
			and simplex.filtration_value <= time
		):
			if _num_colours_in_bitmask(simplex.colours) == 1:
				colour = plot_colours[_bitmask_to_colours(simplex.colours)[0]]
			else:
				colour = "grey"
			ax1.fill(points[0, simplex.vertices], points[1, simplex.vertices], c=colour, alpha=0.2)

	ax1.set_aspect("equal")
	# ax.set_xlabel('Time = ' + f"{0.0:.4f}")

	return ax1


def animate_filtration(
	K: FilteredComplex,
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.float64]],
	filtration_times: Sequence[float],
	animation_length: float,
) -> animation.FuncAnimation:
	"""
	Create animation of 2-skeleton of filtered simplicial complex.

	Args:
		K                : A filtered complex.
		points           : The vertices of ``K`` as a :math:`2\\times N` numpy matrix.
		filtration_times : Sequence of filtration times for which to draw animation frames.
		animation_length : Total length of the animation in seconds.
	"""
	if len(points.shape) != 2:
		raise NotImplementedError

	fig: Figure
	ax: Axes
	fig, ax = plt.subplots()
	plot_colours = np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])
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
			colour = "black"
			alpha = 0.2
		lines[idx] = ax.plot([], [], c=colour, alpha=alpha, linewidth=1)[0]

	for idx, simplex in K.simplices[2].items():
		if _num_colours_in_bitmask(simplex.colours) == 1:
			i = _bitmask_to_colours(simplex.colours)[0]
			colour = plot_colours[i]
		else:
			colour = "grey"
		patches[idx] = ax.fill([], [], c=colour, alpha=0.2)[0]

	ax.set_aspect("equal")
	ax.set_xlabel("Time = " + f"{0.0:.4f}")

	def update(time: float) -> None:
		for idx, simplex in K.simplices[1].items():
			if simplex.filtration_value <= time:
				lines[idx].set_xdata(points[0, simplex.vertices])
				lines[idx].set_ydata(points[1, simplex.vertices])

		for idx, simplex in K.simplices[2].items():
			if simplex.filtration_value <= time:
				patches[idx].set_xy(points[:, simplex.vertices].T)

		ax.set_xlabel(f"Time = {time:.4f}")

	interval = int(np.round(animation_length * 1000 / len(filtration_times)))
	return animation.FuncAnimation(fig=fig, func=update, frames=filtration_times, interval=interval)


def _get_truncation(entrance_times: list[float]) -> float:
	return (
		ub * 2
		if (m := max(entrance_times))
		> 2
		* (
			ub := np.percentile(
				entrance_times,
				[
					95,
				],
			)[0]
		)
		else m
	)

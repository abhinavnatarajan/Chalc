"""Plotting and visualisation utilities."""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from collections.abc import Set as AbstractSet
from itertools import product
from typing import TYPE_CHECKING, Any, Literal, get_args

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import animation
from matplotlib.legend import Legend
from pandas import DataFrame

from .sixpack import (
	DiagramEnsemble,
	SimplexPairings,
	_bitmask_to_colours,
	_colours_are_subset,
	_colours_to_bitmask,
	_num_colours_in_bitmask,
)

if TYPE_CHECKING:
	from matplotlib.artist import Artist
	from matplotlib.axes import Axes
	from matplotlib.figure import Figure

	from chalc.filtration import FilteredComplex

plt.rcParams["animation.html"] = "jshtml"


def plot_sixpack(
	dgms: DiagramEnsemble,
	truncations: Mapping[DiagramEnsemble.DiagramName, float] | float | None = None,
	dimensions: AbstractSet[int] | int | None = None,
	tolerance: float = 0,
) -> tuple[Figure, np.ndarray[tuple[Literal[2], Literal[3]], np.dtype[Axes]]]:
	"""Plot the 6-pack of persistence diagrams returned by :func:`compute <.sixpack.compute>`.

	Args:
		dgms        : The 6-pack of persistence diagrams.
		truncations :
			The maximum entrance time upto which features are plotted.
			Can be either a single value, or a mapping from diagram names to values.
			A sensible default will be calculated (for each individual diagram) if not provided.
		dimensions  :
			The homological dimensions for which to plot features.
			If not provided, all dimensions will be included in the plots.
		tolerance   : Only features with persistence greater than this value will be plotted.

	"""
	diagram_names = get_args(DiagramEnsemble.DiagramName.__value__)
	if dimensions is None:
		if dgms:
			dimensions = set(range(max(dgms.dimensions)))
			# The relative diagram may contain features
			# one dimension higher than that of the filtration.
			dims = {
				name: (dimensions if name != "rel" else dimensions | {max(dimensions) + 1})
				for name in diagram_names
			}
		else:
			# if the diagrams are all empty then we gracefully fail
			# by plotting empty diagrams
			dims = {name: set() for name in diagram_names}
	elif isinstance(dimensions, int):
		dims = {name: {dimensions} for name in diagram_names}
	elif isinstance(dimensions, AbstractSet) and all(isinstance(dim, int) for dim in dimensions):
		dims = {name: dimensions for name in diagram_names}
	else:
		errmsg = "Invalid dimensions argument."
		raise ValueError(errmsg)

	# In general, the dimension of a feature represented by a simplex pair (a, b)
	# is dim(a), unless the diagram is in the kernel, in which case
	# it is dim(a) - 1.
	dim_shift = {name: (1 if name == "ker" else 0) for name in diagram_names}

	if truncations is None:
		truncations = {}
	if isinstance(truncations, Mapping) and all(
		isinstance(val, float) for val in truncations.values()
	):
		truncs = {}
		for diagram_name in diagram_names:
			if diagram_name in truncations and isinstance(truncations[diagram_name], float):
				truncs[diagram_name] = truncations[diagram_name]
				continue
			ycoords = [
				dgms.entrance_times[death_simplex]
				for birth_simplex, death_simplex in dgms[diagram_name].paired
				if dgms.dimensions[birth_simplex] - dim_shift[diagram_name] in dims[diagram_name]
			]
			truncs[diagram_name] = _get_truncation(ycoords) if ycoords else 1.0
	elif isinstance(truncations, float):
		truncs = {diagram_name: truncations for diagram_name in diagram_names}
	else:
		errmsg = "Invalid truncations argument."
		raise TypeError(errmsg)

	plot_pos = {
		name: (val[1], val[0])
		for name, val in zip(diagram_names, product((0, 1, 2), (0, 1)), strict=True)
	}
	plot_titles = {
		"ker": "Kernel",
		"cok": "Cokernel",
		"im": "Image",
		"dom": "Domain",
		"cod": "Codomain",
		"rel": "Relative",
	}
	points_legend = {name: (name == "rel") for name in diagram_names}
	fig, axes = plt.subplots(
		nrows=2,
		ncols=3,
		figsize=[3.5 * 3, 3.5 * 2],
		sharex=False,
		sharey=False,
	)
	for diagram_name in diagram_names:
		_plot_diagram(
			diagram=dgms[diagram_name],
			entrance_times=list(dgms.entrance_times),
			dimensions=list(dgms.dimensions),
			truncation=truncs[diagram_name],
			dims=dims[diagram_name],
			ax=axes[plot_pos[diagram_name]],
			dim_shift=dim_shift[diagram_name],
			title=plot_titles[diagram_name],
			points_legend=points_legend[diagram_name],
			lines_legend=False,
			tolerance=tolerance,
		)
	legends = {}
	for ax in fig.axes:
		ax.tick_params(labelleft=True, labelbottom=True)
		handles, labels = ax.get_legend_handles_labels()
		legends |= dict(zip(labels, handles, strict=True))
	axes[plot_pos["rel"]].get_legend().remove()
	fig.legend(handles=legends.values(), loc="lower center", ncol=5)
	return fig, axes


def plot_diagram(
	dgms: DiagramEnsemble,
	diagram_name: DiagramEnsemble.DiagramName,
	truncation: float | None = None,
	dimensions: AbstractSet[int] | int | None = None,
	ax: Axes | None = None,
	tolerance: float = 0,
) -> Axes | None:
	"""Plot a specific diagram from a 6-pack.

	Args:
		dgms         : The 6-pack of persistence diagrams.
		diagram_name : One of ``'ker'``, ``'cok'``, ``'dom'``, ``'cod'``, ``'im'``, or ``'rel'``.
		truncation   :
			The maximum entrance time for which the diagrams are plotted.
			A sensible default will be calculated if not provided.
		dimensions   :
			The homological dimensions for which to plot features.
			If not provided, all dimensions will be included in the plots.
		ax           :
			A matplotlib axes object.
			If provided then the diagram will be plotted on the given axes.
		tolerance    : Only features with persistence greater than this value will be plotted.

	"""
	dgm = dgms[diagram_name]

	# In general, the dimension of a feature represented by a simplex pair (a, b)
	# is dim(a), unless the diagram is in the kernel, in which case
	# it is dim(a) - 1.
	dim_shift = 1 if diagram_name == "ker" else 0

	if dimensions is None:
		dimensions = set(range(max(dgms.dimensions))) if dgm else set()
	elif isinstance(dimensions, int):
		dimensions = {
			dimensions,
		}

	if truncation is None:
		ycoord = [
			dgms.entrance_times[death_simplex]
			for birth_simplex, death_simplex in dgms[diagram_name].paired
			if dgms.dimensions[birth_simplex] - dim_shift in dimensions
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
		entrance_times=list(dgms.entrance_times),
		dimensions=list(dgms.dimensions),
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
	*,
	diagram: SimplexPairings,
	entrance_times: Sequence[float],  # entrance times of the simplices by index
	dimensions: Sequence[int],  # dimensions of the simplices by index
	truncation: float,  # max birth/death time to plot
	dims: AbstractSet[int],  # dimensions for which features are plotted
	ax: Axes,  # axes to plot on
	title: str,  # title of the plot
	points_legend: bool,  # show legend for points (which dimension)
	lines_legend: bool,  # show legend for lines (truncation, infinity)
	dim_shift: int,  # dimension shift if any; only relevant for kernel
	tolerance: float,
	**kwargs: Any,
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
		for birth_idx, death_idx in diagram.paired
		if dimensions[birth_idx] - dim_shift in dims
		and entrance_times[death_idx] - entrance_times[birth_idx] > tolerance
	] + [
		(
			entrance_times[birth_idx],
			inf_level,
			f"$H_{{{dimensions[birth_idx] - dim_shift}}}$",
		)
		for birth_idx in diagram.unpaired
		if dimensions[birth_idx] - dim_shift in dims
	]
	plot_df = DataFrame.from_records(data=all_pts, columns=["Birth", "Death", "Dimension"])
	plot_df["Dimension"] = plot_df["Dimension"].astype("category")
	ret_ax = sns.scatterplot(
		data=plot_df,
		x="Birth",
		y="Death",
		hue="Dimension",
		ax=ax,
		legend=points_legend,
		**kwargs,
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
		[0, truncate_level],
		[truncate_level, truncate_level],
		"k:",
		alpha=0.4,
		label="Truncated",
	)
	handle.plot([0, inf_level], [inf_level, inf_level], "k--", alpha=0.4, label=r"$\infty$")
	handle.plot([0, 0], [0, inf_level], "k-", alpha=0.4)
	if lines_legend:
		handle.legend()


def draw_filtration(
	K: FilteredComplex,  # noqa: N803
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.floating]],
	time: float,
	include_colours: Collection[int] | None = None,
	ax: Axes | None = None,
) -> Axes:
	r"""Visualise a 2D filtration at given time, optionally including only certain colours.

	Args:
		K               : A filtered 2-dimensional simplicial complex.
		points          : The vertices of ``K`` as a :math:`2\times N` numpy matrix.
		time            : Filtration times for which to draw simplices.
		include_colours :
			Optional collection of colours to include.
			If not specified then all colours will be drawn.
		ax              :
			A matplotlib axes object.
			If provided then the diagram will be plotted on the given axes.

	"""
	if len(points.shape) != 2:  # noqa: PLR2004
		raise NotImplementedError

	if include_colours is None:
		include_colours = {
			_bitmask_to_colours(vertex.colours)[0] for vertex in K.simplices[0].values()
		}

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
	for simplex in K.simplices[1].values():
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

	for simplex in K.simplices[2].values():
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
	K: FilteredComplex,  # noqa: N803
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.floating]],
	filtration_times: Sequence[float],
	animation_length: float,
) -> animation.FuncAnimation:
	r"""Create animation of 2-skeleton of filtered simplicial complex.

	Args:
		K                : A filtered complex.
		points           : The vertices of ``K`` as a :math:`2\times N` numpy matrix.
		filtration_times : Sequence of filtration times for which to draw animation frames.
		animation_length :
			Total length of the animation in seconds, unrelated to the
		filtration times.

	"""
	if len(points.shape) != 2:  # noqa: PLR2004
		raise NotImplementedError

	fig: Figure
	ax: Axes
	fig, ax = plt.subplots()
	plot_colours = np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])
	vertex_colours = []
	for simplex in K.simplices[0].values():
		i = _bitmask_to_colours(simplex.colours)[0]
		vertex_colours.append(i)

	ax.scatter(points[0, :], points[1, :], c=list(plot_colours[vertex_colours]), s=10)

	lines = {}
	patches = {}
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

	def update(time: float) -> list[Artist]:
		modified_artists = []
		for idx, simplex in K.simplices[1].items():
			if simplex.filtration_value <= time:
				lines[idx].set_xdata(points[0, simplex.vertices])
				lines[idx].set_ydata(points[1, simplex.vertices])
				modified_artists.append(lines[idx])

		for idx, simplex in K.simplices[2].items():
			if simplex.filtration_value <= time:
				patches[idx].set_xy(points[:, simplex.vertices].T)
				modified_artists.append(patches[idx])

		ax.set_xlabel(f"Time = {time:.4f}")
		return modified_artists

	interval = int(np.round(animation_length * 1000 / len(filtration_times)))
	return animation.FuncAnimation(
		fig=fig,
		func=update,
		frames=filtration_times,
		interval=interval,
	)  # pyright : ignore


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

"""Plotting and visualisation utilities."""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from itertools import product
from numbers import Real
from typing import TYPE_CHECKING, Any, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from matplotlib.legend import Legend
from pandas import DataFrame

from chalc.sixpack import DiagramName, SimplexPairings, SixPack

if TYPE_CHECKING:
	from matplotlib.artist import Artist
	from matplotlib.axes import Axes
	from matplotlib.figure import Figure

	from chalc.filtration import Filtration

plt.rcParams["animation.html"] = "jshtml"


def plot_sixpack(
	dgms: SixPack,
	truncations: Mapping[DiagramName, float] | float | None = None,
	dimensions: Collection[int] | int | None = None,
	threshold: float = 0,
) -> tuple[Figure, np.ndarray[tuple[Literal[2], Literal[3]], np.dtype[Any]]]:
	"""
	Plot the 6-pack of persistence diagrams returned by :func:`compute <.sixpack.compute>`.

	Args:
		dgms        : The 6-pack of persistence diagrams.
		truncations :
			The maximum entrance time upto which features are plotted.
			Can be either a single value, or a mapping from diagram names to values.
			A sensible default will be calculated (for each individual diagram) if not provided.
		dimensions  :
			The homological dimensions for which to plot features.
			If not provided, all dimensions will be included in the plots.
		threshold   : Only features with persistence greater than this value will be plotted.

	"""
	diagram_names = SixPack.names()
	dims: dict[DiagramName, list[int]]
	if dimensions is None:
		if dgms:
			dim_range = list(range(max(dgms.dimensions.tolist())))
			# The relative diagram may contain features
			# one dimension higher than that of the filtration.
			dims = {
				name: (dim_range if name != "rel" else [*dim_range, max(dim_range) + 1])
				for name in diagram_names
			}
		else:
			# if the diagrams are all empty then we gracefully fail
			# by plotting empty diagrams
			dims = {name: [] for name in diagram_names}
	elif isinstance(dimensions, int):
		dims = {name: [dimensions] for name in diagram_names}
	elif isinstance(dimensions, Collection) and all(isinstance(dim, int) for dim in dimensions):
		dims = {key: sorted(dimensions) for key in diagram_names}
	else:
		errmsg = "Invalid dimensions argument."
		raise ValueError(errmsg)

	# In general, the dimension of a feature represented by a simplex pair (a, b)
	# is dim(a), unless the diagram is in the kernel, in which case
	# it is dim(a) - 1.
	dim_shift: dict[DiagramName, int] = {
		name: (1 if name == "ker" else 0) for name in diagram_names
	}

	if truncations is None:
		truncations = {}
	if isinstance(truncations, Mapping) and all(
		isinstance(val, Real) for val in truncations.values()
	):
		truncs: dict[DiagramName, float] = {}
		for diagram_name in diagram_names:
			if diagram_name in truncations:
				truncs[diagram_name] = truncations[diagram_name]
				continue
			births = []
			deaths = []
			for pair in dgms[diagram_name].paired:
				if (
					dgms.dimensions[pair[0]] - dim_shift[diagram_name] in dims[diagram_name]
					and dgms.entrance_times[pair[1]] - dgms.entrance_times[pair[0]] > threshold
				):
					births.append(dgms.entrance_times[pair[0]].item())
					deaths.append(dgms.entrance_times[pair[1]].item())
			births.extend(
				dgms.entrance_times[unpaired].item()
				for unpaired in dgms[diagram_name].unpaired
				if dgms.dimensions[unpaired] - dim_shift[diagram_name] in dims[diagram_name]
			)
			if not births or not deaths:
				# The diagram is either (no births) or has only infinite features (no deaths).
				truncs[diagram_name] = 1.0
			else:
				truncs[diagram_name] = _get_truncation(births, deaths)
	elif isinstance(truncations, Real):
		truncs = dict.fromkeys(diagram_names, float(truncations))
	else:
		errmsg = "Invalid truncations argument."
		raise TypeError(errmsg)

	plot_pos: dict[DiagramName, tuple[int, int]] = {
		name: (val[1], val[0])
		for name, val in zip(diagram_names, product((0, 1, 2), (0, 1)), strict=True)
	}
	plot_titles: dict[DiagramName, str] = {
		"ker": "Kernel",
		"cok": "Cokernel",
		"im": "Image",
		"dom": "Domain",
		"cod": "Codomain",
		"rel": "Relative",
	}
	points_legend: dict[DiagramName, bool] = {name: (name == "rel") for name in diagram_names}
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
			entrance_times=dgms.entrance_times.tolist(),
			dimensions=dgms.dimensions.tolist(),
			truncation=truncs[diagram_name],
			dims=dims[diagram_name],
			ax=axes[plot_pos[diagram_name]],
			dim_shift=dim_shift[diagram_name],
			title=plot_titles[diagram_name],
			points_legend=points_legend[diagram_name],
			lines_legend=False,
			threshold=threshold,
		)
	legends = {}
	for ax in fig.axes:
		ax.tick_params(labelleft=True, labelbottom=True)
		handles, labels = ax.get_legend_handles_labels()
		legends |= dict(zip(labels, handles, strict=True))
	fig.legend(
		handles=[legends[key] for key in sorted(legends)],
		labels=sorted(legends),
		loc="lower center",
		ncol=5,
	)
	return fig, axes


def plot_diagram(
	dgms: SixPack,
	diagram_name: DiagramName,
	truncation: float | None = None,
	dimensions: Collection[int] | int | None = None,
	ax: Axes | None = None,
	threshold: float = 0,
) -> Axes | None:
	"""
	Plot a specific diagram from a 6-pack.

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
		threshold    : Only features with persistence greater than this value will be plotted.

	"""
	dgm = dgms[diagram_name]

	# In general, the dimension of a feature represented by a simplex pair (a, b)
	# is dim(a), unless the diagram is in the kernel, in which case
	# it is dim(a) - 1.
	dim_shift = 1 if diagram_name == "ker" else 0

	if dimensions is None:
		dimensions = list(range(max(dgms.dimensions.tolist()))) if dgm else []
	elif isinstance(dimensions, int):
		dimensions = [dimensions]
	elif isinstance(dimensions, Collection) and all(isinstance(dim, int) for dim in dimensions):
		dimensions = sorted(dimensions)
	else:
		errmsg = "Invalid dimensions argument."
		raise ValueError(errmsg)

	if truncation is None:
		births = []
		deaths = []
		for pair in dgms[diagram_name].paired:
			if (
				dgms.dimensions[pair[0]] - dim_shift in dimensions
				and dgms.entrance_times[pair[1]] - dgms.entrance_times[pair[0]] > threshold
			):
				births.append(dgms.entrance_times[pair[0]].item())
				deaths.append(dgms.entrance_times[pair[1]].item())

		births.extend(
			dgms.entrance_times[unpaired].item()
			for unpaired in dgms[diagram_name].unpaired
			if dgms.dimensions[unpaired] - dim_shift in dimensions
		)
		truncation = 1.0 if not births or not deaths else _get_truncation(births, deaths)

	if ax is None:
		ax1: Axes
		_, ax1 = plt.subplots()
	else:
		ax1 = ax

	titles: dict[DiagramName, str] = {
		"ker": "Kernel",
		"cok": "Cokernel",
		"im": "Image",
		"dom": "Domain",
		"cod": "Codomain",
		"rel": "Relative",
	}

	_plot_diagram(
		diagram=dgm,
		entrance_times=dgms.entrance_times.tolist(),
		dimensions=dgms.dimensions.tolist(),
		truncation=truncation,
		dims=dimensions,
		ax=ax1,
		title=titles[diagram_name],
		points_legend=True,
		lines_legend=True,
		dim_shift=dim_shift,
		threshold=threshold,
	)
	return ax1


def _plot_diagram(
	*,
	diagram: SimplexPairings,
	entrance_times: Sequence[float],  # entrance times of the simplices by index
	dimensions: Sequence[int],  # dimensions of the simplices by index
	truncation: float,  # max birth/death time to plot
	dims: Collection[int],  # dimensions for which features are plotted
	ax: Axes,  # axes to plot on
	title: str,  # title of the plot
	points_legend: bool,  # show legend for points (which dimension)
	lines_legend: bool,  # show legend for lines (truncation, infinity)
	dim_shift: int,  # dimension shift if any; only relevant for kernel
	threshold: float,
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
			feature_dimension,
		)
		for birth_idx, death_idx in diagram.paired
		if (feature_dimension := dimensions[birth_idx] - dim_shift) in dims
		and entrance_times[death_idx] - entrance_times[birth_idx] > threshold
	] + [
		(
			entrance_times[birth_idx],
			inf_level,
			feature_dimension,
		)
		for birth_idx in diagram.unpaired
		if (feature_dimension := dimensions[birth_idx] - dim_shift) in dims
	]
	plot_colours = np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])
	plot_df = DataFrame.from_records(data=all_pts, columns=["Birth", "Death", "Dimension"])
	plot_df["Dimension"] = plot_df["Dimension"].astype("category")
	for d in sorted(plot_df["Dimension"].cat.categories):
		ax.scatter(
			plot_df.loc[plot_df["Dimension"] == d, "Birth"],
			plot_df.loc[plot_df["Dimension"] == d, "Death"],
			c=plot_colours[d],
			edgecolors="white",
			label=f"$H_{{{d}}}$",
			**kwargs,
		)
	ax.set(xlabel=None)
	ax.set(ylabel=None)
	ax.set(title=title)
	ax.set_aspect("equal")
	has_legend = any(True for c in ax.get_children() if isinstance(c, Legend))
	if points_legend and has_legend:
		ax.legend(loc="lower right")
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
	K: Filtration,  # noqa: N803
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.floating]],
	time: float,
	include_colours: Collection[int] | None = None,
	ax: Axes | None = None,
	plot_colours: str | np.ndarray | None = None,
) -> Axes:
	r"""
	Visualise a 2D filtration at given time, optionally including only certain colours.

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
		plot_colours    :
			Either the name of a matplotlib qualitative colour map, or a list of colours.
			If not provided, the current matplotlib colour cycle will be used.

	"""
	if len(points.shape) != 2:  # noqa: PLR2004
		raise NotImplementedError

	all_colours = {vertex.colours[0] for vertex in K.simplices[0].values()}
	colours_map = {old: new for new, old in enumerate(sorted(all_colours))}

	include_colours = all_colours if include_colours is None else include_colours

	ax1 = plt.subplots()[1] if ax is None else ax

	if isinstance(plot_colours, str):
		num_colours = K.dimension - 1
		colour_indices = np.arange(num_colours)
		plot_colours = plt.get_cmap(plot_colours)(colour_indices)
	elif plot_colours is None:
		plot_colours = np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])

	# Plot the vertices
	vertex_plot_colours = []
	vertices_to_plot = []
	for idx in range(K.num_vertices):
		vertex = K.simplices[0][idx]
		if (set(vertex.colours).issubset(include_colours)) and vertex.filtration_value <= time:
			vertex_plot_colours.append(plot_colours[colours_map[vertex.colours[0]]])
			vertices_to_plot.append(idx)

	ax1.scatter(
		points[0, vertices_to_plot],
		points[1, vertices_to_plot],
		c=list(vertex_plot_colours),
		s=10,
	)

	# Plot the edges
	for edge in K.simplices[1].values():
		if set(edge.colours).issubset(include_colours) and edge.filtration_value <= time:
			colour, alpha = (
				(plot_colours[colours_map[edge.colours[0]]], 0.5)
				if len(edge.colours) == 1
				else ("black", 0.2)
			)
			ax1.plot(
				points[0, edge.vertices],
				points[1, edge.vertices],
				c=colour,
				alpha=alpha,
				linewidth=1,
			)

	for triangle in K.simplices[2].values():
		if set(triangle.colours).issubset(include_colours) and triangle.filtration_value <= time:
			colour = (
				plot_colours[colours_map[triangle.colours[0]]]
				if len(triangle.colours) == 1
				else "grey"
			)
			ax1.fill(
				points[0, triangle.vertices],
				points[1, triangle.vertices],
				c=colour,
				alpha=0.2,
			)

	ax1.set_aspect("equal")

	return ax1


def animate_filtration(
	K: Filtration,  # noqa: N803
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.floating]],
	filtration_times: Sequence[float],
	animation_length: float,
	plot_colours: str | np.ndarray | None = None,
) -> animation.FuncAnimation:
	r"""
	Create animation of 2-skeleton of filtered simplicial complex.

	Args:
		K                : A filtered complex.
		points           : The vertices of ``K`` as a :math:`2\times N` numpy matrix.
		filtration_times : Sequence of filtration times for which to draw animation frames.
		animation_length :
			Total length of the animation in seconds, unrelated to the
			filtration times.
		plot_colours    :
			Either the name of a matplotlib qualitative colour map, or a list of colours.
			If not provided, the current matplotlib colour cycle will be used.

	"""
	if len(points.shape) != 2:  # noqa: PLR2004
		raise NotImplementedError

	fig: Figure
	ax: Axes
	fig, ax = plt.subplots()

	all_colours = {vertex.colours[0] for vertex in K.simplices[0].values()}
	colours_map = {old: new for new, old in enumerate(sorted(all_colours))}

	if isinstance(plot_colours, str):
		num_colours = K.dimension - 1
		colour_indices = np.arange(num_colours)
		plot_colours = plt.get_cmap(plot_colours)(colour_indices)
	elif plot_colours is None:
		plot_colours = np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])

	vertex_plot_colours = [
		plot_colours[colours_map[K.simplices[0][i].colours[0]]] for i in range(K.num_vertices)
	]

	ax.scatter(points[0, :], points[1, :], c=vertex_plot_colours, s=10)

	lines = {}
	patches = {}
	for idx, edge in K.simplices[1].items():
		colour, alpha = (
			(plot_colours[colours_map[edge.colours[0]]], 0.5)
			if len(edge.colours) == 1
			else ("black", 0.2)
		)
		lines[idx] = ax.plot([], [], c=colour, alpha=alpha, linewidth=1)[0]

	for idx, triangle in K.simplices[2].items():
		colour = (
			plot_colours[colours_map[triangle.colours[0]]] if len(triangle.colours) == 1 else "grey"
		)
		patches[idx] = ax.fill([], [], c=colour, alpha=0.2)[0]

	ax.set_aspect("equal")
	ax.set_xlabel("Time = " + f"{0.0:.4f}")

	def update(time: float) -> list[Artist]:
		modified_artists = []
		for idx, edge in K.simplices[1].items():
			if edge.filtration_value <= time:
				lines[idx].set_xdata(points[0, edge.vertices])
				lines[idx].set_ydata(points[1, edge.vertices])
				modified_artists.append(lines[idx])

		for idx, triangle in K.simplices[2].items():
			if triangle.filtration_value <= time:
				patches[idx].set_xy(points[:, triangle.vertices].T)
				modified_artists.append(patches[idx])

		ax.set_xlabel(f"Time = {time:.4f}")
		return modified_artists

	interval = int(np.round(animation_length * 1000 / len(filtration_times)))
	return animation.FuncAnimation(
		fig=fig,
		func=update,
		frames=filtration_times,
		interval=interval,
	)


def _get_truncation(births: list[float], deaths: list[float]) -> float:
	perc95_death = np.percentile(deaths, [95])[0]
	max_death = max(deaths)
	death_sensible_trunc = perc95_death * 2 if max_death > 2 * perc95_death else max_death
	return max(max(births) * 1.05, death_sensible_trunc)

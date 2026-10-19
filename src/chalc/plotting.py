"""Plotting and visualisation utilities."""

from __future__ import annotations

from collections.abc import Callable, Collection, Mapping, Sequence
from numbers import Integral
from typing import TYPE_CHECKING, Any, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

if TYPE_CHECKING:
	from matplotlib.artist import Artist
	from matplotlib.axes import Axes
	from matplotlib.figure import Figure

	from chalc.filtration import Filtration, Simplex
	from chalc.morphisms import (
		AugmentedPersistenceModule,
		IndexArray,
		PersistenceModuleName,
		SixPack,
	)

plt.rcParams["animation.html"] = "jshtml"

__all__ = ["animate_filtration", "draw_filtration", "plot_persistence_diagram", "plot_sixpack"]

type _1DArrayFloat = Sequence[float] | np.ndarray[tuple[int], np.dtype[np.floating]]  # noqa: PYI042


def plot_sixpack(
	sixpack: SixPack,
	truncations: Mapping[PersistenceModuleName, float | None] | float | None = None,
	degrees: Collection[int] | IndexArray | int | np.integer | None = None,
	threshold: float = 0,
) -> tuple[Figure, np.ndarray[tuple[Literal[2], Literal[3]], np.dtype[Any]]]:
	"""
	Plot the persistence diagrams for a :meth:`chalc.morphisms.FiltrationMorphism.sixpack`.

	Args:
		sixpack        : The 6-pack of persistence modules.
		truncations :
			The maximum entrance time up to which features are plotted.
			Can be either a single value, or a mapping from persistence module
			names to values. A sensible default is calculated for module names
			omitted from the mapping.
		degrees  :
			The homological degrees for which to plot features.
			If not provided, all degrees will be included in the plots.
		threshold   : Only features with persistence greater than this value will be plotted.

	Returns:
		The figure and a ``(2, 3)`` array of axes. The top row contains domain,
		codomain, and image; the bottom row contains kernel, cokernel, and relative.
	"""
	if isinstance(truncations, Mapping):
		truncations = dict(truncations)
	else:
		truncations = dict.fromkeys(sixpack.module_names(), truncations)

	fig, axes = plt.subplots(
		nrows=2,
		ncols=3,
		figsize=[3.5 * 3, 3.5 * 2],
		sharex=False,
		sharey=False,
		constrained_layout=True,
	)

	plot_pos: dict[PersistenceModuleName, tuple[int, int]] = {
		"domain": (0, 0),
		"kernel": (1, 0),
		"codomain": (0, 1),
		"cokernel": (1, 1),
		"image": (0, 2),
		"relative": (1, 2),
	}
	for name, module in sixpack:
		plot_persistence_diagram(
			module,
			truncation=truncations.get(name),
			degrees=degrees,
			threshold=threshold,
			points_legend=False,
			lines_legend=False,
			ax=axes[plot_pos[name]],
		)
		axes[plot_pos[name]].set_title(name.title())
	legends = {}
	for ax in fig.axes:
		ax.tick_params(labelleft=True, labelbottom=True)
		handles, labels = ax.get_legend_handles_labels()
		legends |= dict(zip(labels, handles, strict=True))
	if legends:
		fig.legend(
			handles=[legends[key] for key in sorted(legends)],
			labels=sorted(legends),
			loc="lower center",
			ncol=5,
		)
	return fig, axes


def plot_persistence_diagram(
	persistence_module: AugmentedPersistenceModule,
	*,
	truncation: float | None = None,
	degrees: Collection[int] | IndexArray | int | np.integer | None = None,
	ax: Axes | None = None,
	threshold: float = 0,
	kwargs_per_dim: Mapping[int, dict] | None = None,
	points_legend: bool = True,
	lines_legend: bool = True,
	**kwargs,
) -> Axes:
	"""
	Plot the persistence diagram for a persistence module.

	Args:
		persistence_module:
			The augmented persistence module.
		truncation:
			The maximum entrance time for which the persistence module is plotted.
			A positive finite default is calculated from the retained intervals.
		degrees:
			The homological degrees in which to plot features.
			If not provided, all degrees will be included in the plots.
		ax:
			A matplotlib axes object.
			If provided then the persistence diagram will be plotted on the given axes.
		threshold:
			Only features with persistence greater than this value will be plotted.
			The threshold is applied before truncating their displayed coordinates.
		kwargs_per_dim:
			Keyword arguments for each degree, as a mapping `{ degree : dictionary of kwargs }`.
			These will be passed to the matplotlib scatter function, and take precedence over
			any other keyword arguments.
		points_legend:
			Whether to show a legend for the feature degrees.
		lines_legend:
			Whether to show a legend for the truncation lines.
		**kwargs:
			Additional keyword arguments to pass to the matplotlib scatter function.

	Returns:
		The axes containing the diagram, creating a new figure if ``ax`` is omitted.
	"""
	degrees = _select_degrees(persistence_module, degrees)

	intervals_by_degree = {
		d: mat[(mat[:, 1] - mat[:, 0]) > threshold]
		for d, mat in persistence_module.intervals(degrees).items()
	}
	if truncation is None:
		matrices = list(intervals_by_degree.values())
		intervals = np.concatenate(matrices) if matrices else np.empty((0, 2))
		truncation = _get_truncation(intervals[:, 0], intervals[:, 1])
	if not np.isfinite(truncation) or truncation <= 0:
		msg = "Truncation must be a positive finite number."
		raise ValueError(msg)

	if ax is None:
		_, ax = plt.subplots()
	truncate_level = truncation * 1.04
	inf_level = truncation * 1.08

	def truncate(values: np.ndarray) -> np.ndarray:
		return np.where(
			np.isinf(values),
			inf_level,
			np.where(values > truncation, truncate_level, values),
		)

	marker_colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
	for d, mat in intervals_by_degree.items():
		if mat.size == 0:
			continue
		style = _scatter_options(
			marker_colours[d % len(marker_colours)],
			kwargs,
			kwargs_per_dim.get(d, {}) if kwargs_per_dim is not None else {},
		)
		ax.scatter(truncate(mat[:, 0]), truncate(mat[:, 1]), label=f"$H_{{{d}}}$", **style)
	ax.set(xlabel=None, ylabel=None)
	ax.set_aspect("equal")
	ax.plot([0, inf_level], [0, inf_level], "k-", alpha=0.4)
	ax.plot(
		[0, truncate_level],
		[truncate_level, truncate_level],
		"k:",
		alpha=0.4,
		label="Truncated" if lines_legend else "_nolegend_",
	)
	ax.plot(
		[0, inf_level],
		[inf_level, inf_level],
		"k--",
		alpha=0.4,
		label=r"$\infty$" if lines_legend else "_nolegend_",
	)
	ax.plot([0, 0], [0, inf_level], "k-", alpha=0.4)
	handles, labels = ax.get_legend_handles_labels()
	selected = [
		(h, label)
		for h, label in zip(handles, labels, strict=True)
		if points_legend or not label.startswith("$H_{")
	]
	if selected:
		legend_handles, legend_labels = zip(*selected, strict=True)
		ax.legend(legend_handles, legend_labels)
	return ax


def _select_degrees(
	module: AugmentedPersistenceModule,
	degrees: Collection[int] | IndexArray | int | np.integer | None,
) -> IndexArray:
	if degrees is None:
		return module.nonempty_degrees()
	if isinstance(degrees, Integral):
		degrees = [int(degrees)]
	if not isinstance(degrees, Collection) or not all(isinstance(d, Integral) for d in degrees):
		msg = "Invalid degrees argument."
		raise ValueError(msg)
	return np.intersect1d(module.nonempty_degrees(), np.fromiter(degrees, dtype=np.int64))


def _scatter_options(colour: str, global_options: dict, degree_options: dict) -> dict:
	style = {"color": colour, "edgecolors": "white", **global_options}
	if "c" in style:
		style.pop("color", None)
	if "c" in degree_options:
		style.pop("color", None)
	elif "color" in degree_options:
		style.pop("c", None)
	return style | degree_options


def draw_filtration(
	filtration: Filtration,
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.floating]],
	select: Callable[[Simplex, int], bool] | None = None,
	*,
	ax: Axes | None = None,
	plot_colours: str | np.ndarray | None = None,
	vertices_kwargs: dict | None = None,
	edges_kwargs: dict | None = None,
	faces_kwargs: dict | None = None,
) -> Axes:
	r"""
	Draw selected vertices, edges, and triangles of a filtration in the plane.

	Args:
		filtration      : A filtered simplicial complex. Only its 2-skeleton is drawn.
		points          : The vertices of ``filtration`` as a :math:`2\times N` numpy matrix.
		select          :
			Optional callback ``select(simplex, index)`` returning ``True`` if the
			simplex should be drawn. It receives a simplex handle and its zero-based
			index in ``filtration.boundary_matrix()``. The callback is called once
			for each simplex of dimension at most 2, in boundary-matrix order.
			Simplices of higher dimension still count towards the indices. If omitted,
			all vertices, edges, and triangles are drawn.
		ax              :
			A matplotlib axes object.
			If provided then the filtration will be plotted on the given axes.
		plot_colours    :
			Either the name of a matplotlib qualitative colour map, or a list of colours.
			If not provided, the current matplotlib colour cycle will be used.
		vertices_kwargs :
			Keyword arguments to pass to the matplotlib scatter function for plotting points.
		edges_kwargs    :
			Keyword arguments to pass to the matplotlib plot function for plotting edges.
		faces_kwargs    :
			Keyword arguments to pass to the matplotlib fill function for plotting faces.

	Returns:
		The axes on which the simplices were drawn.

	Raises:
		ValueError: If ``points`` does not have shape ``(2, filtration.num_vertices)``.

	Notes:
		Selecting a simplex does not automatically select its faces or vertices.
		The callback must not modify ``filtration``.
		Representative supports must use the indices of ``filtration.boundary_matrix()``.
		For representatives computed from a dimension-limited boundary matrix,
		pass the corresponding skeleton as ``filtration``.
		Representatives from a quotient use mapping-cylinder indices and must
		first be mapped to the plotted filtration's indices.

	"""
	if points.shape != (2, filtration.num_vertices):
		msg = "points must have shape (2, filtration.num_vertices)."
		raise ValueError(msg)

	simplex_handles = filtration.simplices
	all_colours = sorted({vertex.colours[0] for vertex in simplex_handles[0].values()})
	colours_map = {colour: index for index, colour in enumerate(all_colours)}

	if isinstance(plot_colours, str):
		plot_colours = plt.get_cmap(plot_colours)(np.arange(len(all_colours)))
	elif plot_colours is None:
		plot_colours = np.asarray(plt.rcParams["axes.prop_cycle"].by_key()["color"])

	ax1 = plt.subplots()[1] if ax is None else ax
	vertex_plot_colours = []
	vertices_to_plot = []
	style: dict[str, Any]
	for index, (boundary, label, _, colours) in enumerate(filtration.boundary_matrix()):
		dimension = len(boundary) - 1 if boundary else 0
		simplex = simplex_handles[dimension][label]
		if dimension > 2 or (select is not None and not select(simplex, index)):  # noqa: PLR2004
			continue

		vertices = simplex.vertices
		monochromatic = len(colours) == 1
		colour = plot_colours[colours_map[colours[0]]] if monochromatic else None
		if dimension == 0:
			vertices_to_plot.append(vertices[0])
			vertex_plot_colours.append(colour)
		elif dimension == 1:
			style = {
				"linewidth": 1,
				"color": colour if monochromatic else "black",
				"alpha": 0.5 if monochromatic else 0.3,
				**(edges_kwargs or {}),
			}
			ax1.plot(
				points[0, vertices],
				points[1, vertices],
				**style,
			)
		else:
			style = {
				"edgecolor": None,
				"facecolor": colour if monochromatic else "grey",
				"alpha": 0.2,
				**(faces_kwargs or {}),
			}
			ax1.fill(
				points[0, vertices],
				points[1, vertices],
				**style,
			)

	if vertices_to_plot:
		style = {"s": 10, "c": vertex_plot_colours, **(vertices_kwargs or {})}
		if "color" in style:
			style.pop("c", None)
		ax1.scatter(
			points[0, vertices_to_plot],
			points[1, vertices_to_plot],
			**style,
		)

	ax1.set_aspect("equal")

	return ax1


def animate_filtration(
	filtration: Filtration,
	points: np.ndarray[tuple[Literal[2], int], np.dtype[np.floating]],
	filtration_times: Sequence[float],
	*,
	animation_length: float,
	plot_colours: str | np.ndarray | None = None,
	vertices_kwargs: dict | None = None,
	edges_kwargs: dict | None = None,
	faces_kwargs: dict | None = None,
) -> animation.FuncAnimation:
	r"""
	Create animation of 2-skeleton of filtered simplicial complex.

	Args:
		filtration       : A filtered complex.
		points           : The vertices of ``filtration`` as a :math:`2\times N` numpy matrix.
		filtration_times : Sequence of filtration times for which to draw animation frames.
		animation_length :
			Total length of the animation in seconds, unrelated to the
			filtration times.
		plot_colours     :
			Either the name of a matplotlib qualitative colour map, or a list of colours.
			If not provided, the current matplotlib colour cycle will be used.
		vertices_kwargs  :
			Keyword arguments to pass to the matplotlib scatter function for plotting points.
		edges_kwargs     :
			Keyword arguments to pass to the matplotlib plot function for plotting edges.
		faces_kwargs     :
			Keyword arguments to pass to the matplotlib fill function for plotting faces.

	"""
	if len(points.shape) != 2:  # noqa: PLR2004
		raise NotImplementedError
	if vertices_kwargs is None:
		vertices_kwargs = {}
	if edges_kwargs is None:
		edges_kwargs = {}
	if faces_kwargs is None:
		faces_kwargs = {}

	fig: Figure
	ax: Axes
	fig, ax = plt.subplots()

	all_colours = {vertex.colours[0] for vertex in filtration.simplices[0].values()}
	colours_map = {old: new for new, old in enumerate(sorted(all_colours))}

	if isinstance(plot_colours, str):
		num_colours = filtration.dimension - 1
		colour_indices = np.arange(num_colours)
		plot_colours = plt.get_cmap(plot_colours)(colour_indices)
	if plot_colours is None:
		plot_colours = np.array(plt.rcParams["axes.prop_cycle"].by_key()["color"])

	vertex_plot_colours = [
		plot_colours[colours_map[filtration.simplices[0][i].colours[0]]]
		for i in range(filtration.num_vertices)
	]

	vertices_kwargs = {"s": 10, **vertices_kwargs}
	ax.scatter(
		points[0, :],
		points[1, :],
		c=list(vertex_plot_colours),
		**vertices_kwargs,
	)

	edges_kwargs = {"linewidth": 1, **edges_kwargs}
	lines = {}
	for idx, edge in filtration.simplices[1].items():
		if len(edge.colours) == 1:
			colour = plot_colours[colours_map[edge.colours[0]]]
			alpha = 0.5
		else:
			colour = "black"
			alpha = 0.3
		lines[idx] = ax.plot(
			[],
			[],
			color=colour,
			alpha=alpha,
			**edges_kwargs,
		)[0]

	faces_kwargs = {"alpha": 0.2, **faces_kwargs}
	patches = {}
	for idx, triangle in filtration.simplices[2].items():
		colour = (
			plot_colours[colours_map[triangle.colours[0]]] if len(triangle.colours) == 1 else "grey"
		)
		patches[idx] = ax.fill(
			[],
			[],
			edgecolor=None,
			facecolor=colour,
			**faces_kwargs,
		)[0]

	ax.set_aspect("equal")
	ax.set_xlabel("Time = " + f"{0.0:.4f}")

	def update(time: float) -> list[Artist]:
		modified_artists: list[Artist] = []
		for idx, edge in filtration.simplices[1].items():
			if edge.filtration_value <= time:
				lines[idx].set_xdata(points[0, edge.vertices])
				lines[idx].set_ydata(points[1, edge.vertices])
				modified_artists.append(lines[idx])

		for idx, triangle in filtration.simplices[2].items():
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


def _get_truncation(births: _1DArrayFloat, deaths: _1DArrayFloat) -> float:
	finite_births = np.asarray(births)[np.isfinite(births)]
	finite_deaths = np.asarray(deaths)[np.isfinite(deaths)]
	birth_bound = float(np.max(finite_births)) * 1.05 if finite_births.size else 0.0
	death_bound = 0.0
	if finite_deaths.size:
		perc95 = float(np.percentile(finite_deaths, 95))
		death_bound = min(float(np.max(finite_deaths)), 2 * perc95)
	return max(0.0, birth_bound, death_bound) or 1.0

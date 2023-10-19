from __future__ import annotations

from matplotlib.figure import Figure
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

def plot_sixpack(
	dgms                  : DiagramEnsemble,
	*,
	truncation            : float | None = None,
	max_diagram_dimension : int | None   = None
	) -> tuple[Figure, Annotated[NDArray, "Axes"]] :
	
	"""
	Plots the 6-pack of persistence diagrams returned by :func:`compute <.sixpack.compute>`.

	Args:
		dgms : The 6-pack of persistence diagrams.

	Keyword Args:
		truncation : The maximum entrance time for which the diagrams are plotted. A sensible default will be calculated if not provided.
		max_diagram_dimension : The maximum homological dimension for which to plot points. If not provided, all dimensions will be included in the plots.
	"""
	if truncation is None:
		truncation = _get_truncation(dgms.entrance_times)
	if max_diagram_dimension is None:
		max_diagram_dimension = max(dgms.dimensions)
	fig, axes = plt.subplots(nrows=2, ncols=3, figsize=[3.5 * 3, 3.5 * 2], sharex = True, sharey = True)
	def applied_plot(
		dgm,
		max_dim = max_diagram_dimension,
		**kwargs
		):
		_plot_diagram(
			dgm,
			dgms.entrance_times,
			dgms.dimensions,
			truncation,
			max_dim = max_dim,
			**kwargs
		)

	applied_plot(
		dgms.ker,
		ax        = axes[0, 0],
		title     = "Kernel",
		dim_shift = 1,
	)
	applied_plot(
		dgms.rel,
		ax            = axes[0, 1],
		title         = "Relative",
		max_dim       = max_diagram_dimension + 1,
		points_legend = True
	)
	applied_plot(
		dgms.cok,
		ax    = axes[0, 2],
		title = "Cokernel",
	)
	applied_plot(
		dgms.dom,
		ax    = axes[1, 0],
		title = "Domain",
	)
	applied_plot(
		dgms.im,
		ax    = axes[1, 1],
		title = "Image",
	)
	applied_plot(
		dgms.cod,
		ax    = axes[1, 2],
		title = "Codomain",
	)
	legends = dict()
	for ax in fig.axes:
		handles, labels = ax.get_legend_handles_labels()
		for h,l in zip(handles, labels):
			legends[l] = h
	axes[0, 1].get_legend().remove()
	fig.legend(handles=legends.values(), loc = 'upper center', bbox_to_anchor = (0.5, 0.05), ncol = 5)
	return fig, axes

@interpolate_docstring()
def plot_diagram(
	dgms                  : DiagramEnsemble,
	diagram_name          : str,
	*,
	truncation            : float | None = None,
	max_diagram_dimension : int | None   = None ,
	) -> tuple[Figure, Annotated[NDArray, "Axes"]] :
	"""
	Plot a specific diagram from a 6-pack.

	Args:
		dgms : The 6-pack of persistence diagrams.
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
		ax            = ax,
		dim_shift     = 1 if diagram_name == 'ker' else 0,
		points_legend = True,
		lines_legend  = True,
		title         = titles[diagram_name],
		max_dim       = max_diagram_dimension
	)
	return fig, ax

def _plot_diagram(
	diagram,
	entrance_times,
	dimensions,
	truncation,
	max_dim       = 2,
	ax            = None,
	title         = None,
	points_legend = False,
	lines_legend  = False,
	dim_shift     = 0,
	**kwargs
	) -> None :
	# Truncate all times
	truncate_level = truncation * 1.04
	inf_level = truncation * 1.08
	entrance_times = [
		et if et <= truncation else truncate_level if et is not np.inf else inf_level for et in entrance_times
	]
	all_pts = [
		(
			entrance_times[birth_idx],
			entrance_times[death_idx],
			f"$H_{{{dimensions[birth_idx] - dim_shift}}}$",
		)
		for birth_idx, death_idx in diagram.paired
		if entrance_times[death_idx] != entrance_times[birth_idx]
		and dimensions[birth_idx] - dim_shift <= max_dim
	] + \
	[
		(
			entrance_times[birth_idx],
			inf_level,
			f"$H_{{{dimensions[birth_idx] - dim_shift}}}$",
		)
		for birth_idx in diagram.unpaired
		if dimensions[birth_idx] - dim_shift <= max_dim
	]
	df = pd.DataFrame(data=all_pts, columns=["Birth", "Death", "Dimension"])
	df["Dimension"] = df["Dimension"].astype("category")
	ret_ax = sns.scatterplot(
		data   = df,
		x      = "Birth",
		y      = "Death",
		hue    = "Dimension",
		ax     = ax,
		legend = points_legend,
		**kwargs
	)
	ret_ax.set(xlabel = None)
	ret_ax.set(ylabel = None)
	ret_ax.set(title  = title)
	ret_ax.set_aspect('equal')
	if points_legend:
		sns.move_legend(ret_ax, "lower right")
	handle = ax if ax is not None else plt
	handle.plot([0, inf_level], [0, inf_level], "k-", alpha=0.4)
	handle.plot([0, inf_level], [truncate_level, truncate_level], 'k:', alpha = 0.4, label = "Truncated")
	handle.plot([0, inf_level], [inf_level, inf_level], 'k--', alpha = 0.4, label = "$\\infty$")
	handle.plot([0, 0], [0, inf_level], "k-", alpha=0.4)
	if lines_legend:
		handle.legend()

def draw_filtration(
	K               : FilteredComplex,
	points          : NDArray[np.float64],
	time            : float,
	*,
	include_colours : list[int] | None = None
	) -> tuple[Figure, Axes] :
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

	fig : Figure
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
	K                : FilteredComplex,
	points           : NDArray[np.float64],
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

	fig : Figure
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
	return ub  if (m := max(entrance_times)) > 1.5 * (ub := scipy.stats.scoreatpercentile(entrance_times, 95)) else m

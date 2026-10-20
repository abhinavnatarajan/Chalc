"""Regression tests for plotting the revised persistence API."""

from collections.abc import Iterator

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import PathCollection
from matplotlib.colors import to_rgba
from matplotlib.patches import Polygon

from chalc import chromatic
from chalc.filtration import Filtration, Simplex, complete_complex
from chalc.morphisms import (
	AugmentedPersistenceModule,
	KChromaticGluing,
	SixPack,
	SixPackPairings,
	SubChromaticInclusion,
)
from chalc.plotting import draw_filtration, plot_persistence_diagram, plot_sixpack


@pytest.fixture(autouse=True)
def close_figures() -> Iterator[None]:
	yield
	plt.close("all")


def test_plot_diagram_degrees_and_styles() -> None:
	module = AugmentedPersistenceModule({0: 2, 1: None}, [0.0, 0.0, 2.0], {0: 0, 1: 12})
	styles = {12: {"c": "red", "s": 40}}
	_, ax = plt.subplots()
	assert plot_persistence_diagram(module, ax=ax, kwargs_per_dim=styles, edgecolors="black") is ax
	assert len(ax.collections) == 2
	assert ax.get_legend() is not None
	assert styles == {12: {"c": "red", "s": 40}}
	assert np.isfinite(ax.collections[1].get_offsets()).all()
	assert (
		# pyrefly: ignore [bad-argument-type]
		plot_persistence_diagram(module, degrees=np.int64(0)).collections[0].get_label()
		== "$H_{0}$"
	)
	assert len(plot_persistence_diagram(module, degrees={12}).collections) == 1
	assert len(plot_persistence_diagram(module, degrees=[1]).collections) == 0


@pytest.mark.parametrize(
	"module",
	[AugmentedPersistenceModule(), AugmentedPersistenceModule({0: None}, [0.0], {0: 0})],
)
def test_plot_empty_or_infinite(module: AugmentedPersistenceModule) -> None:
	ax = plot_persistence_diagram(module)
	ax.figure.canvas.draw()
	assert np.isfinite(ax.get_xlim()).all()
	assert np.isfinite(ax.get_ylim()).all()


def test_plot_threshold_before_truncation() -> None:
	module = AugmentedPersistenceModule(
		{0: 3, 1: None, 2: 2},
		[0.0, 3.0, 4.0, 10.0],
		{0: 0, 1: 0, 2: 0},
	)
	ax = plot_persistence_diagram(module, truncation=1.0, threshold=2.0)
	# Both long-lived intervals survive even though displayed persistence is small.
	np.testing.assert_allclose(
		np.asarray(ax.collections[0].get_offsets()),
		[[0.0, 1.04], [1.04, 1.08]],
	)
	assert (
		len(
			np.asarray(
				plot_persistence_diagram(module, threshold=10.0).collections[0].get_offsets(),
			),
		)
		== 1
	)


@pytest.mark.parametrize("points", [False, True])
@pytest.mark.parametrize("lines", [False, True])
def test_plot_legend_flags(*, points: bool, lines: bool) -> None:
	module = AugmentedPersistenceModule({0: None}, [0.0], {0: 0})
	ax = plot_persistence_diagram(module, points_legend=points, lines_legend=lines)
	legend = ax.get_legend()
	labels = {text.get_text() for text in legend.get_texts()} if legend else set()
	assert ("$H_{0}$" in labels) == points
	assert ("Truncated" in labels) == lines
	assert (r"$\infty$" in labels) == lines


@pytest.mark.parametrize("empty", [False, True])
def test_plot_sixpack(*, empty: bool) -> None:
	pack = (
		SixPack()
		if empty
		else SixPack(pairings=SixPackPairings(domain={0: None}), entrance_times=[0.0], degrees=[0])
	)
	fig, axes = plot_sixpack(pack, truncations={"domain": 2.0})
	assert axes.shape == (2, 3)
	assert [[ax.get_title() for ax in row] for row in axes] == [
		["Domain", "Codomain", "Image"],
		["Kernel", "Cokernel", "Relative"],
	]
	fig.canvas.draw()
	assert all(np.isfinite(ax.get_xlim()).all() for ax in fig.axes)
	if not empty:
		assert [text.get_text() for text in fig.legends[0].get_texts()] == ["$H_{0}$"]


def test_plot_module_with_negative_times() -> None:
	module = AugmentedPersistenceModule({0: 1}, [-2.0, -1.0], {0: 0})
	ax = plot_persistence_diagram(module)
	ax.figure.canvas.draw()
	np.testing.assert_array_equal(ax.collections[0].get_offsets(), [[-2.0, -1.0]])


def test_plotting_accepts_degree_arrays() -> None:
	module = AugmentedPersistenceModule({0: 2, 1: None}, [0.0, 0.0, 2.0], {0: 0, 1: 12})
	ax = plot_persistence_diagram(module, degrees=np.array([0, 12], dtype=np.int64))
	assert len(ax.collections) == 2


@pytest.fixture
def drawing_filtration() -> tuple[Filtration, np.ndarray]:
	filtration = Filtration(5, 3)
	filtration.add_simplex([0, 1, 2, 3], 1.0)
	filtration.add_simplex([3, 4], 2.0)
	for label, vertex in filtration.simplices[0].items():
		vertex.set_colour(label % 2)
	filtration.propagate_colours()
	points = np.array([[0.0, 1.0, 0.0, 1.0, 2.0], [0.0, 0.0, 1.0, 1.0, 0.0]])
	return filtration, points


@pytest.mark.parametrize("include_faces", [False, True])
def test_draw_callback_receives_full_boundary_indices(
	drawing_filtration: tuple[Filtration, np.ndarray],
	*,
	include_faces: bool,
) -> None:
	filtration, points = drawing_filtration
	matrix = filtration.boundary_matrix()
	seen = []

	def select(simplex: Simplex, index: int) -> bool:
		boundary, label, time, colours = matrix[index]
		assert simplex.dimension == max(0, len(boundary) - 1)
		assert simplex.label == label
		assert simplex.filtration_value == time
		assert simplex.colours == colours
		assert simplex is filtration.simplices[simplex.dimension][label]
		seen.append(index)
		return index == len(matrix) - 1

	_, ax = plt.subplots()
	assert draw_filtration(filtration, points, select, include_faces=include_faces, ax=ax) is ax
	assert seen == [i for i, column in enumerate(matrix) if len(column[0]) <= 3]
	assert len(matrix[-2][0]) == 4  # The tetrahedron precedes the selected edge.
	assert len(ax.lines) == 1
	np.testing.assert_array_equal(ax.lines[0].get_xydata(), points[:, [3, 4]].T)
	if include_faces:
		assert len(ax.collections) == 1
		np.testing.assert_array_equal(ax.collections[0].get_offsets(), points[:, [3, 4]].T)
	else:
		assert not ax.collections
	assert not ax.patches


@pytest.mark.parametrize("dimension", [0, 1, 2, 3])
@pytest.mark.parametrize("include_faces", [False, True])
def test_draw_default_selection_and_low_dimensions(dimension: int, *, include_faces: bool) -> None:
	filtration = complete_complex(4, dimension)
	points = np.array([[0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, 1.0]])
	ax = draw_filtration(filtration, points, include_faces=include_faces, plot_colours="tab10")
	assert len(ax.collections) == 1
	np.testing.assert_array_equal(ax.collections[0].get_offsets(), points.T)
	assert len(ax.lines) == (6 if dimension >= 1 else 0)
	assert len(ax.patches) == (4 if dimension >= 2 else 0)


@pytest.mark.parametrize("include_faces", [False, True])
def test_draw_empty_selection(
	drawing_filtration: tuple[Filtration, np.ndarray],
	*,
	include_faces: bool,
) -> None:
	filtration, points = drawing_filtration
	ax = draw_filtration(
		filtration, points, lambda _simplex, _index: False, include_faces=include_faces
	)
	assert not ax.collections
	assert not ax.lines
	assert not ax.patches


@pytest.mark.parametrize("include_faces", [False, True])
def test_draw_selected_triangle_and_its_faces(*, include_faces: bool) -> None:
	filtration = Filtration(6, 2)
	filtration.add_simplex([0, 1, 2], 1.0)
	# A later, distant triangle catches reuse of the last simplex's geometry.
	filtration.add_simplex([3, 4, 5], 2.0)
	points = np.array([[0.0, 1.0, 0.0, 10.0, 0.0, 10.0], [0.0, 0.0, 1.0, 0.0, 10.0, 10.0]])
	ax = draw_filtration(
		filtration,
		points,
		lambda simplex, _: simplex.vertices == [0, 1, 2],
		include_faces=include_faces,
	)
	assert len(ax.patches) == 1
	triangle = ax.patches[0]
	assert isinstance(triangle, Polygon)
	np.testing.assert_array_equal(triangle.get_xy(), points[:, [0, 1, 2, 0]].T)
	if include_faces:
		assert len(ax.lines) == 3
		assert {frozenset(map(tuple, np.asarray(line.get_xydata()))) for line in ax.lines} == {
			frozenset(map(tuple, points[:, vertices].T)) for vertices in ([0, 1], [0, 2], [1, 2])
		}
		assert len(ax.collections) == 1
		np.testing.assert_array_equal(ax.collections[0].get_offsets(), points[:, [0, 1, 2]].T)
	else:
		assert not ax.lines
		assert not ax.collections


def test_draw_time_and_colour_selector(drawing_filtration: tuple[Filtration, np.ndarray]) -> None:
	filtration, points = drawing_filtration
	ax = draw_filtration(
		filtration,
		points,
		lambda simplex, _: simplex.filtration_value <= 1.0 and set(simplex.colours) <= {0},
	)
	np.testing.assert_array_equal(ax.collections[0].get_offsets(), points[:, [0, 2, 4]].T)
	assert len(ax.lines) == 1
	np.testing.assert_array_equal(ax.lines[0].get_xydata(), points[:, [0, 2]].T)
	assert not ax.patches


def test_draw_style_overrides(drawing_filtration: tuple[Filtration, np.ndarray]) -> None:
	filtration, points = drawing_filtration
	vertices = {"color": "purple", "s": 30}
	edges = {"color": "crimson", "alpha": 1.0, "linewidth": 3}
	faces = {"facecolor": "lime", "alpha": 0.6}
	ax = draw_filtration(
		filtration,
		points,
		vertices_kwargs=vertices,
		edges_kwargs=edges,
		faces_kwargs=faces,
	)
	collection = ax.collections[0]
	assert isinstance(collection, PathCollection)
	np.testing.assert_allclose(np.asarray(collection.get_facecolor()), [to_rgba("purple")])
	assert collection.get_sizes().tolist() == [30]
	assert all(line.get_color() == "crimson" for line in ax.lines)
	assert all(line.get_alpha() == 1.0 and line.get_linewidth() == 3 for line in ax.lines)
	assert all(patch.get_facecolor() == to_rgba("lime", 0.6) for patch in ax.patches)
	assert vertices == {"color": "purple", "s": 30}
	assert edges == {"color": "crimson", "alpha": 1.0, "linewidth": 3}
	assert faces == {"facecolor": "lime", "alpha": 0.6}


@pytest.mark.parametrize("shape", [(3, 5), (2, 4), (5,), (2, 5, 1)])
def test_draw_requires_one_2d_point_per_vertex(
	drawing_filtration: tuple[Filtration, np.ndarray],
	shape: tuple[int, ...],
) -> None:
	filtration, _ = drawing_filtration
	with pytest.raises(ValueError, match="points must have shape"):
		# pyrefly: ignore [bad-argument-type]
		draw_filtration(filtration, np.zeros(shape))


@pytest.mark.parametrize("gluing", [False, True])
def test_draw_kernel_representatives(*, gluing: bool) -> None:
	points = np.array([[0.0, 1.0, 1.0, 0.0, 0.5], [0.0, 0.0, 1.0, 1.0, 0.5]])
	filtration = chromatic.alpha(points, [0, 0, 0, 0, 1], max_num_threads=1)
	morphism = KChromaticGluing(filtration, 1) if gluing else SubChromaticInclusion(filtration, {0})
	pack = morphism.sixpack(get_representatives=True)
	assert pack.representatives is not None
	matrix = filtration.boundary_matrix()
	if gluing:
		assert pack.metadata is not None
		index_map = pack.metadata["domain_index_to_codomain_index"]
	else:
		index_map = range(len(matrix))
	cycles = [
		{index_map[int(index)] for index in cycle}
		for cycle in pack.representatives.kernel.values()
		if len(matrix[index_map[int(cycle[0])]][0]) == 2
	]
	assert len(cycles) == 1
	support = cycles[0]
	assert len(support) == 4
	boundary = set()
	for index in support:
		boundary ^= set(matrix[index][0])
	assert not boundary

	ax = draw_filtration(
		filtration,
		points,
		lambda _, index: index in support,
		edges_kwargs={"color": "crimson", "alpha": 1.0},
	)
	expected_edges = {frozenset(edge) for edge in ((0, 1), (1, 2), (2, 3), (3, 0))}
	assert {
		frozenset(filtration.simplices[1][matrix[index][1]].vertices) for index in support
	} == expected_edges
	assert len(ax.lines) == 4
	for line, index in zip(ax.lines, sorted(support), strict=True):
		vertices = filtration.simplices[1][matrix[index][1]].vertices
		np.testing.assert_array_equal(line.get_xydata(), points[:, vertices].T)
	assert len(ax.collections) == 1
	np.testing.assert_array_equal(ax.collections[0].get_offsets(), points[:, [0, 1, 2, 3]].T)
	assert not ax.patches

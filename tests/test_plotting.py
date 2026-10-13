"""Regression tests for plotting the revised persistence API."""

from collections.abc import Iterator

import matplotlib.pyplot as plt
import numpy as np
import pytest

from chalc.morphisms import AugmentedPersistenceModule, SixPack
from chalc.plotting import plot_persistence_diagram, plot_sixpack


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
	np.testing.assert_allclose(ax.collections[0].get_offsets(), [[0.0, 1.04], [1.04, 1.08]])
	assert len(plot_persistence_diagram(module, threshold=10.0).collections[0].get_offsets()) == 1


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
	pack = SixPack() if empty else SixPack(domain={0: None}, entrance_times=[0.0], degrees=[0])
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

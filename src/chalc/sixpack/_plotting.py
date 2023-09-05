from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from numbers import Real
from ._types import DiagramEnsemble

def plot_sixpack(dgms : DiagramEnsemble,
                 truncation : Real, 
                 max_diagram_dim : int = 2) -> None :
    """
    Plots the 6-pack of persistence diagrams returned by :func:`compute`.

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
        plot_diagram(
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
    return fig


def plot_diagram(
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
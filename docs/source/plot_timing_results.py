import json
import logging
from pathlib import Path

import chromatic_tda as chro
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import chalc as ch

funcs = {
	"chromatic alpha multi-threaded": ch.chromatic.alpha,
	"chromatic alpha single-threaded": lambda points, colours: ch.chromatic.alpha(
		points,
		colours,
		max_num_threads=1,
	),
	"chromatic delcech multi-threaded": ch.chromatic.delcech,
	"chromatic delcech single-threaded": lambda points, colours: ch.chromatic.delcech(
		points,
		colours,
		max_num_threads=1,
	),
	"chromatic delrips multi-threaded": ch.chromatic.delrips,
	"chromatic delrips single-threaded": lambda points, colours: ch.chromatic.delrips(
		points,
		colours,
		max_num_threads=1,
	),
	"chromatic-tda": lambda points, colours: chro.ChromaticAlphaComplex(points.T, colours),
}

logger = logging.getLogger()

# Load the data from the JSON file
output_path = Path("/home/abhinav/Documents/software/chalc/docs/source/timing_results.json")
with Path(output_path).open("r") as f:
	experiment = json.load(f)
	logger.info("Loaded experiment result from json.")

# Create plots of the runtimes vs number of points, for each function and number of colours
df = pd.DataFrame(experiment)

unique_num_colours = sorted(df["num_colours"].unique())

# Get unique function families and assign colors from matplotlib's default cycler
unique_function_families = []
for name in funcs:
	if "single-threaded" in name:
		base_name = name.replace(" single-threaded", "")
	elif "multi-threaded" in name:
		base_name = name.replace(" multi-threaded", "")
	else:
		base_name = name
	if base_name not in unique_function_families:
		unique_function_families.append(base_name)

# Create color mapping using matplotlib's default color cycle
plt.style.use("seaborn-v0_8-darkgrid")

# Set font sizes for all plot elements
plt.rcParams.update({
    "font.size": 16,           # Default font size
    "axes.titlesize": 18,      # Title font size
    "axes.labelsize": 16,      # Axis label font size
    "xtick.labelsize": 14,     # X-axis tick label font size
    "ytick.labelsize": 14,     # Y-axis tick label font size
    "legend.fontsize": 16,     # Legend font size
    "figure.titlesize": 18     # Figure title font size
})

default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
function_colors = {
	family: default_colors[i % len(default_colors)]
	for i, family in enumerate(unique_function_families)
}

# Create separate figures for each number of colours
fig, axes = plt.subplots(len(unique_num_colours), 1, figsize=(12, 12))
if len(axes) == 1:
	axes = [axes]
for i, num_colours in enumerate(unique_num_colours):
	ax = axes[i]
	colour_data = df[df["num_colours"] == num_colours]
	# Plot each function
	for name in np.unique(colour_data["name"]):
		func_data = colour_data[colour_data["name"] == name]
		# Determine base function name and line style
		if "single-threaded" in name:
			base_name = name.replace(" single-threaded", "")
			linestyle = ":"
		elif "multi-threaded" in name:
			base_name = name.replace(" multi-threaded", "")
			linestyle = "-"
		else:
			base_name = name
			linestyle = "-"
		# Get color for this function family
		color = function_colors.get(base_name, "black")
		ax.plot(
			func_data["num_points"],
			func_data["mean_exec_time"],
			marker="o",
			label=name,
			linewidth=2,
			markersize=6,
			color=color,
			linestyle=linestyle,
		)
		ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", handlelength=3)
	ax.set_xlabel("Number of Points")
	ax.set_ylabel("Mean Execution Time (s)")
	ax.set_title(f"Runtime vs Number of Points ({num_colours} colours)")
	ax.set_yscale("log")

plt.tight_layout()

# Save each plot separately
if __file__:
	plot_path = Path(__file__).parent / "timing_results.png"
	plt.savefig(plot_path, dpi=300, bbox_inches="tight")
	logmsg = f"Plot saved to {plot_path}"
	logger.info(logmsg)

plt.show()


# Also create a summary table
summary_table = df.pivot_table(
	values="mean_exec_time", index=["num_points", "num_colours"], columns="name", aggfunc="mean"
)

print("\nSummary Table (Mean Execution Time in seconds):")  # noqa: T201
print(summary_table.to_string(float_format=lambda x: f"{x:.4f}"))  # noqa: T201

plt.show()

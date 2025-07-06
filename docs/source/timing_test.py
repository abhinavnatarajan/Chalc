import functools
import itertools
import json
import logging
import sys
import timeit
from pathlib import Path

import chromatic_tda as chro
import numpy as np

import chalc as ch

funcs = {
	"chromatic alpha multi-threaded": lambda points, colours: ch.chromatic.alpha(
		points,
		colours,
		max_num_threads=8,
	),
	"chromatic alpha single-threaded": lambda points, colours: ch.chromatic.alpha(
		points,
		colours,
		max_num_threads=1,
	),
	"chromatic delcech multi-threaded": lambda points, colours: ch.chromatic.delcech(
		points,
		colours,
		max_num_threads=8,
	),
	"chromatic delcech single-threaded": lambda points, colours: ch.chromatic.delcech(
		points,
		colours,
		max_num_threads=1,
	),
	"chromatic delrips multi-threaded": lambda points, colours: ch.chromatic.delrips(
		points,
		colours,
		max_num_threads=8,
	),
	"chromatic delrips single-threaded": lambda points, colours: ch.chromatic.delrips(
		points,
		colours,
		max_num_threads=1,
	),
	"chromatic-tda": lambda points, colours: chro.ChromaticAlphaComplex(points.T, colours),
}

experiment = []
num_trials = 3
num_points_list = [100, 250, 500, 1000, 2500, 5000]
num_colours_list = [2, 3]

logger = logging.getLogger()
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

for num_points, num_colours in itertools.product(num_points_list, num_colours_list):
	rng = np.random.default_rng(40)
	points = rng.uniform(size=(2, num_points))
	colours = rng.integers(0, num_colours, size=num_points).tolist()
	for name, func in funcs.items():
		apply_func = functools.partial(func, points, colours)
		timer = timeit.Timer(apply_func)
		logmsg = f"Timing {name} with {num_points} points and {num_colours} colours..."
		logger.info(logmsg)
		times = timer.repeat(repeat=num_trials, number=1)
		logmsg = f"mean={np.mean(times):.4f}s, var={np.var(times)}s"
		logger.info(logmsg)
		experiment.append(
			{
				"name": name,
				"num_points": num_points,
				"num_colours": len(set(colours)),
				"mean_exec_time": np.mean(times),
				"var_exec_time": np.var(times),
			},
		)

# Save experiment results to JSON file
if __file__:
	output_path = Path(__file__).parent / "timing_results.json"
	with Path(output_path).open("w") as f:
		json.dump(experiment, f, indent=2)
	logmsg = f"Results saved to {output_path}."
	logger.info(logmsg)

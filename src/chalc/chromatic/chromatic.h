/*
This file is part of Chalc.

Chalc: Chromatic Alpha Complexes.
Based on: di Montesano et. al., “Persistent Homology of Chromatic Alpha
Complexes”. Online preprint available at http://arxiv.org/abs/2212.03128.
Accessed: 2023-02-28 22:07:23 UTC.
DOI: 10.48550/arXiv.2212.03128.

Project homepage: http://abhinavnatarajan.github.io/Chalc

Copyright (c) 2023 Abhinav Natarajan

Contributors:
Abhinav Natarajan

Licensing:
Chalc is released under the GNU General Public License ("GPL").

GNU General Public License ("GPL") copyright permissions statement:
**************************************************************************
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#ifndef CHROMATIC_H
	#define CHROMATIC_H

	#include <CGAL/tags.h>
	#include <Eigen/Dense>
	#include <chalc/filtration/filtration.h>

namespace chalc {

// Compute a Delaunay triangulation from a collection of coordinate vectors
template <typename Concurrency_tag = CGAL::Sequential_tag>
auto delaunay(const Eigen::MatrixXd& X, const std::vector<colour_t>& colours) -> Filtration;

// Compute the chromatic delrips complex
auto delrips(const Eigen::MatrixXd& points, const std::vector<colour_t>& colours) -> Filtration;

// Compute the chromatic delrips complex with parallelisation
auto delrips_parallel(
	const Eigen::MatrixXd&       points,
	const std::vector<colour_t>& colours,
	const int                    max_num_threads
) -> Filtration;

// Compute the chromatic alpha complex
auto alpha(const Eigen::MatrixXd& points, const std::vector<colour_t>& colours)
	-> Filtration;

// Compute the chromatic alpha complex with parallelisation
auto alpha_parallel(
	const Eigen::MatrixXd&       points,
	const std::vector<colour_t>& colours,
	const int                    max_num_threads
) -> Filtration;

// Compute the chromatic Delaunay--Cech complex
auto delcech(const Eigen::MatrixXd& points, const std::vector<colour_t>& colours)
	-> Filtration;

// Compute the chromatic Delaunay--Cech complex with parallelisation
auto delcech_parallel(
	const Eigen::MatrixXd&       points,
	const std::vector<colour_t>& colours,
	const int                    max_num_threads
) -> Filtration;
}  // namespace chalc

#endif

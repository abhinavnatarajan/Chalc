/*
    This file is part of Chalc.

    Chalc: Chromatic Alpha Complexes.
    Based on: di Montesano et. al., “Persistent Homology of Chromatic Alpha
    Complexes”. Online preprint available at http://arxiv.org/abs/2212.03128.
    Accessed: 2023-02-28 22:07:23 UTC.
    DOI: 10.48550/arXiv.2212.03128.

    Project homepage:    http://abhinavnatarajan.github.io/Chalc

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

	#include <Eigen/Dense>
	#include <chalc/common.h>
	#include <chalc/filtration/filtration.h>

namespace chalc {
namespace chromatic {

// Create a Delaunay triangulation from a collection of coordinate vectors
FilteredComplex delaunay(const Eigen::MatrixXd& X, const std::vector<index_t>& colours);

// Create the chromatic delrips complex
FilteredComplex delrips(const Eigen::MatrixXd&      points,
                                          const std::vector<index_t>& colours);

// Create the chromatic alpha complex
std::tuple<FilteredComplex, bool> alpha(const Eigen::MatrixXd&      points,
                                        const std::vector<index_t>& colours);

// Create the chromatic Del-Cech complex
std::tuple<FilteredComplex, bool> delcech(const Eigen::MatrixXd&      points,
                                          const std::vector<index_t>& colours);
}  // namespace chromatic
}  // namespace chalc

#endif

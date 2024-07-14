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
#ifndef CHALC_COMMON_H
	#define CHALC_COMMON_H

	#include <map>
	#include <memory>
	#include <numeric>
	#include <optional>
	#include <tuple>
	#include <vector>

namespace chalc {
typedef double        value_t;
typedef long long int index_t;

namespace stl {
using std::vector, std::map, std::tuple, std::iota, std::tie, std::shared_ptr, std::make_shared,
	std::enable_shared_from_this, std::weak_ptr, std::optional;
}
}  // namespace chalc

#endif

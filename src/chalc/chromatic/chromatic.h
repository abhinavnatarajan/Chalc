#pragma once
#ifndef CHROMATIC_H
	#define CHROMATIC_H

	#include <CGAL/tags.h>
	#include <Eigen/Dense>
	#include <chalc/filtration/filtration.h>

namespace chalc {

// Compute a Delaunay triangulation from a collection of coordinate vectors
template <typename Concurrency_tag = CGAL::Sequential_tag>
auto delaunay(const Eigen::MatrixXd& X, const std::vector<Colour>& colours) -> Filtration;

// Compute the chromatic delrips complex
auto delaunay_rips(const Eigen::MatrixXd& points, const std::vector<Colour>& colours) -> Filtration;

// Compute the chromatic delrips complex with parallelisation
auto delaunay_rips_parallel(
	const Eigen::MatrixXd&       points,
	const std::vector<Colour>& colours,
	const int                    max_num_threads
) -> Filtration;

// Compute the chromatic alpha complex
auto alpha(const Eigen::MatrixXd& points, const std::vector<Colour>& colours)
	-> Filtration;

// Compute the chromatic alpha complex with parallelisation
auto alpha_parallel(
	const Eigen::MatrixXd&       points,
	const std::vector<Colour>& colours,
	const int                    max_num_threads
) -> Filtration;

// Compute the chromatic Delaunay--Čech complex
auto delaunay_cech(const Eigen::MatrixXd& points, const std::vector<Colour>& colours)
	-> Filtration;

// Compute the chromatic Delaunay--Čech complex with parallelisation
auto delaunay_cech_parallel(
	const Eigen::MatrixXd&       points,
	const std::vector<Colour>& colours,
	const int                    max_num_threads
) -> Filtration;
}  // namespace chalc

#endif

#pragma once
#ifndef CHROMATIC_H
#define CHROMATIC_H

#include <Eigen/dense>
#include "filtration.h"

namespace chalc {
    using Eigen::Matrix, Eigen::MatrixXd, Eigen::Matrix2Xd, Eigen::Dynamic;
    using Colouring = vector<size_t>;

    // Stratify a coloured point set
    // Points are provided as columns of a matrix
    // Colours are provided as a vector
    MatrixXd stratify(const MatrixXd& points, const Colouring& colours);

    // Create a Delaunay triangulation from a collection of coordinate vectors
    FilteredComplex delaunay_complex(const MatrixXd& X);

    // Create the weak chromatic alpha complex
    FilteredComplex weak_chromatic_alpha_complex(const MatrixXd& points, const Colouring& colours);

    // Create the chromatic alpha complex
    // FilteredComplex chromatic_alpha_complex(const Matrix<double, 2, Dynamic>& points, const Colouring& colours);

}

#endif
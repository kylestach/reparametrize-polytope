#include "path.h"
#include <iostream>
#include <chrono>

HermiteSplinePath::HermiteSplinePath(WorldSpace p0, WorldSpace d0, WorldSpace p1, WorldSpace d1) {
    Eigen::Matrix<double, 3, 4> target;
    target.block<3, 1>(0, 0) = p0;
    target.block<3, 1>(0, 1) = d0;
    target.block<3, 1>(0, 2) = p1;
    target.block<3, 1>(0, 3) = d1;
    Eigen::Matrix<double, 4, 4> basis;
    basis << 1, 0, -3, 2,
             0, 1, -2, 1,
             0, 0, 3, -2,
             0, 0, -1, 1;
    coefficients = target * basis;
}

PathSample HermiteSplinePath::calculate(double t) {
    PathSample result;

    Eigen::Matrix<double, 4, 1> ts;
    ts << 1, t, t * t, t * t * t;
    result.position = coefficients * ts;

    Eigen::Matrix<double, 4, 1> ts_deriv;
    ts_deriv << 0, 1, 2 * t, 3 * t * t;
    result.direction = coefficients * ts_deriv;

    Eigen::Matrix<double, 4, 1> ts_deriv2;
    ts_deriv2 << 0, 0, 2, 6 * t;
    result.acceleration = coefficients * ts_deriv2;

    if (result.direction.lpNorm<Eigen::Infinity>() > 1e-6) {
        result.curvature = (
                result.direction(0) * result.acceleration(1) -
                result.direction(1) * result.acceleration(0)
                ) / std::pow(
                    result.direction(0) * result.direction(0) + result.direction(1) * result.direction(1)
                    , 1.5);
    } else {
        Eigen::Matrix<double, 4, 1> ts_deriv3;
        ts_deriv3 << 0, 0, 0, 6;
        WorldSpace jerk = coefficients * ts_deriv3;
        result.curvature = (
                result.acceleration(0) * jerk(1) -
                result.acceleration(1) * jerk(0)
                ) / std::pow(
                    result.acceleration(0) * result.acceleration(0) +
                    result.acceleration(1) * result.acceleration(1)
                    , 1.5);
    }

    return result;
}

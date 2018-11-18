#pragma once

#include <iostream>
#include <limits>
#include <tuple>
#include "eigen3/Eigen/Core"
#include "path.h"

// Find the maximum value \ell \in \R such that:
// lb <= A(lt + n) <= ub
template<int NumStates, int NumConstraints>
double maximize_linear(
        const Eigen::Matrix<double, NumConstraints, NumStates>& A,
        const Eigen::Matrix<double, NumStates, 1>& t,
        const Eigen::Matrix<double, NumStates, 1>& n,
        const Eigen::Matrix<double, NumConstraints, 1>& lb,
        const Eigen::Matrix<double, NumConstraints, 1>& ub) {
    auto lb_transformed = lb - A * n;
    auto ub_transformed = ub - A * n;

    double ell_upper = std::numeric_limits<double>::infinity();
    auto At = A * t;
    for (int constraint = 0; constraint < NumConstraints; constraint++) {
        double coordinate = At(constraint, 0);
        double ell_lb = lb_transformed(constraint) / coordinate;
        double ell_ub = ub_transformed(constraint) / coordinate;

        if (ell_lb < 0 && ell_ub < 0) {
            continue;
        }

        if (std::isnan(ell_lb)) {
            ell_lb = -std::numeric_limits<double>::infinity();
        }
        if (std::isnan(ell_ub)) {
            ell_ub = std::numeric_limits<double>::infinity();
        }

        ell_upper = std::min(ell_upper, std::max(ell_lb, ell_ub));
    }

    return ell_upper;
}

WorldSpace maximize_acceleration(
        WorldSpace velocity_initial, double curvature,
        Eigen::Matrix<double, 4, 3> world_to_wheel,
        Eigen::Matrix<double, 4, 1> a_max) {
    double speed_initial = velocity_initial.block<2, 1>(0, 0).norm();
    WorldSpace forward_vector = velocity_initial;
    WorldSpace left_vector(-forward_vector(1), forward_vector(0), 0.0);

    WorldSpace normal_acceleration = curvature * speed_initial * left_vector;

    double ell = maximize_linear<3, 4>(world_to_wheel, forward_vector, normal_acceleration, -a_max, a_max);

    return ell * forward_vector + normal_acceleration;
}

std::tuple<WorldSpace, WorldSpace> maximize_velocity(
        WorldSpace direction, double curvature,
        Eigen::Matrix<double, 4, 3> world_to_wheel,
        Eigen::Matrix<double, 4, 1> a_max) {
    WorldSpace normal_vector(-curvature * direction(1), curvature * direction(0), 0);
    // Find the maximal normal force we can provide
    double ell = maximize_linear<3, 4>(world_to_wheel, normal_vector, WorldSpace::Zero(), -a_max, a_max);
    double normal_acceleration = (normal_vector * ell).norm();
    double speed =  std::sqrt(std::abs(normal_acceleration / curvature));
    WorldSpace max_velocity = direction / direction.norm() * speed;
    return std::make_tuple(max_velocity, normal_vector * ell);
}

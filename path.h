#pragma once
#include "eigen3/Eigen/Core"

using WorldSpace = Eigen::Matrix<double, 3, 1>;

struct PathSample {
    WorldSpace position;
    WorldSpace direction;
    WorldSpace acceleration;
    double curvature;
};

class HermiteSplinePath {
    public:
        HermiteSplinePath(WorldSpace p0, WorldSpace d0, WorldSpace p1, WorldSpace d1);
        PathSample calculate(double t);

    private:
        Eigen::Matrix<double, 3, 4> coefficients;
};

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

#include "reparametrize.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include "eigen3/Eigen/Core"

int main() {
    // Convert world space to wheel space.
    //
    //  1/------\4
    //   |  y   |
    //   |  |_x |
    //   |      |
    //  2\------/3
    //
    // Wheel angles, measured from the x-axis.
    double theta_1 = M_PI * 0.25,
           theta_2 = M_PI * 0.75,
           theta_3 = M_PI * 0.25,
           theta_4 = M_PI * 0.75;
    double r = 0.07;
    Eigen::Matrix<double, 4, 3> M;
    M << std::cos(theta_1), std::sin(theta_1), r * std::sin(theta_1 - M_PI * 0.75),
         std::cos(theta_2), std::sin(theta_2), r * std::sin(theta_2 + M_PI * 0.75),
         std::cos(theta_3), std::sin(theta_3), r * std::sin(theta_3 + M_PI * 0.25),
         std::cos(theta_4), std::sin(theta_4), r * std::sin(theta_4 - M_PI * 0.25);

    // Maximum accelerations in wheel space.
    double a_max = 3.0;

    std::vector<WorldSpace> space(5000);

    double total_time = 0.0;
    int num_iterations = 1;

    std::ofstream file("test.csv");

    WorldSpace velocity(0.0, 0.0, 0.0);

    PathPlanner planner(M, a_max);
    planner.plan(WorldSpace(0, 0, 0), WorldSpace(1, 0, 0), WorldSpace(1, 1, 0), WorldSpace(1, 0, 0));
}

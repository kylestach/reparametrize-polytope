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
    Eigen::Matrix<double, 4, 1> a_max(0.1, 0.1, 0.1, 0.1);

    std::vector<WorldSpace> space(5000);

    double total_time = 0.0;
    int num_iterations = 1;

    std::ofstream file("test.csv");

    double curvature = 1;
    double dt = 0.01;
    WorldSpace position = WorldSpace::Zero();
    for (int j = 0; j < num_iterations; j++) {
        auto start = std::chrono::high_resolution_clock::now();
        space[0] = WorldSpace(.01, 0, 0);
        for (int i = 1; i < space.size(); i++) {
            auto accel = maximize_acceleration(space[i - 1], curvature, M, a_max);

            WorldSpace max_velocity, max_normal_acceleration;
            std::tie(max_velocity, max_normal_acceleration) = maximize_velocity(space[i - 1], curvature, M, a_max);

            auto accelerated_velocity = space[i - 1] + dt * accel;
            auto max_accelerated_velocity = max_velocity + dt * max_normal_acceleration;

            WorldSpace active_acceleration;
            if (max_velocity.norm() > accelerated_velocity.norm()) {
                space[i] = accelerated_velocity;
                active_acceleration = accel;
            } else {
                space[i] = max_accelerated_velocity;
                active_acceleration = max_normal_acceleration;
            }
            if (std::abs((M * active_acceleration)(0)) > 1.01 || std::abs((M * active_acceleration)(1)) > 1.01) {
                std::cout << "FAIL" << std::endl;
            }
            position += (space[i - 1] + space[i]) / 2.0 * dt;
            file << position(0) << "\t" << position(1) << std::endl;
        }
        auto diff = std::chrono::high_resolution_clock::now() - start;
        total_time += std::chrono::duration_cast<std::chrono::milliseconds>(diff).count() / 1000.0;
    }
    file.close();
    std::cout << "Average time: " << total_time / num_iterations << std::endl;
}

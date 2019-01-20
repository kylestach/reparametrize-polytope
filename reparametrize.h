#pragma once

#include <iostream>
#include <limits>
#include <tuple>
#include <vector>
#include "eigen3/Eigen/Core"
#include "path.h"
#include "optimize.hpp"

class PathPlanner {
    public:
        PathPlanner(Eigen::Matrix<double, 4, 3> world_to_wheel,
                    double a_max, double v_max = std::numeric_limits<double>::infinity())
            : world_to_wheel(world_to_wheel),
              a_max(a_max, a_max, a_max, a_max),
              v_max(v_max) {}

        void plan(WorldSpace p0, WorldSpace v0, WorldSpace p1, WorldSpace v1) {
            generate_samples(PathSample{p0, v0, WorldSpace::Zero(), 0.0}, PathSample{p1, v1, WorldSpace::Zero(), 0.0}, 200);
            velocity_pass();
            backwards_pass(v1);
            forwards_pass(v0);
            double t = 0;
            for (size_t i = 0; i < samples.size(); i++) {
                std::cout << t << " " << velocities[i].norm() << " " << samples[i].position(0) << " " << samples[i].position(1) << " " << samples[i].position(2) << std::endl;
                if (i < samples.size() - 1) {
                    t += time_deltas[i];
                }
            }
        }

    private:
        void generate_samples(PathSample start, PathSample end, int num_segments) {
            double directional_length = (start.position - end.position).norm();
            HermiteSplinePath path(
                    start.position,
                    start.direction.normalized() * directional_length,
                    end.position,
                    end.direction.normalized() * directional_length);

            double ds = 1.0 / num_segments;
            for (int i = 0; i < num_segments + 1; i++) {
                double s = i * ds;
                samples.push_back(path.calculate(s));
                velocities.push_back(WorldSpace(0, 0, 0));
            }

            for (int i = 0; i < num_segments; i++) {
                time_deltas.push_back(0);
            }
        }

        void forwards_pass(WorldSpace v0) {
            assert(samples.size() >= 1);

            size_t num_segments = samples.size() - 1;
            velocities[0] = v0;
            for (size_t i = 0; i < num_segments; i++) {
                // Pass over segment i, from sample i to i-1
                PathSample sample_0 = samples[i];
                PathSample sample_1 = samples[i + 1];
                double curvature = (sample_0.curvature + sample_1.curvature) / 2.0;
                WorldSpace velocity_0 = velocities[i];

                Eigen::Matrix<double, 3, 1> direction_0_normalized = sample_0.direction;

                {
                    double direction_0_norm = direction_0_normalized.block<2, 1>(0, 0).norm();
                    if (direction_0_norm > 1e-6) {
                        direction_0_normalized /= direction_0_norm;
                    } else {
                        direction_0_normalized = sample_0.acceleration / sample_0.acceleration.block<2, 1>(0, 0).norm();
                        if (i == num_segments) {
                            direction_0_normalized = -direction_0_normalized;
                        }
                    }
                }

                // Optimize for the maximum forward acceleration, or maintain the maximum
                // velocity that we can keep without having to provide too much centripetal force.
                WorldSpace normal_vector(
                        -std::copysign(1, curvature) * direction_0_normalized(1),
                        std::copysign(1, curvature) * direction_0_normalized(0), 0);
                WorldSpace centripetal_acceleration =
                    normal_vector * velocity_0.dot(velocity_0) * curvature;

                Eigen::Matrix<double, 4, 1> wheel_acceleration;
                double ell;
                bool exact = OptimizeOnLine(
                        Eigen::Matrix<double, 4, 1>(world_to_wheel * direction_0_normalized),
                        Eigen::Matrix<double, 4, 1>(world_to_wheel * centripetal_acceleration),
                        Eigen::Matrix<double, 4, 1>(-a_max),
                        Eigen::Matrix<double, 4, 1>(a_max),
                        &wheel_acceleration,
                        &ell,
                        false);

                if (exact) {
                    WorldSpace forward_acceleration = ell * direction_0_normalized;
                    WorldSpace net_acceleration = centripetal_acceleration + forward_acceleration;

                    // Approximate the segment as a linear piece, with the forward acceleration and velocity components:
                    // dx = a*dt^2 + v*dt
                    // 0 = a*dt^2 + v*dt - dx
                    // dt = (-v +/- sqrt(v^2 + 4*a*dx))/(2a)
                    // dt must be positive, so:
                    // dt = (-v + sqrt(v^2 + 4*a*dx))/(2a)
                    double a_mag = forward_acceleration.norm();
                    double v_mag = velocity_0.norm();
                    double dx = (sample_0.position - sample_1.position).norm();
                    double dt = (-v_mag + std::sqrt(v_mag * v_mag + 4 * a_mag * dx)) / (2 * a_mag);

                    if (dt > time_deltas[i]) {
                        time_deltas[i] = dt;
                        velocities[i + 1] = velocity_0 + dt * net_acceleration;
                    }
                }
            }
        }

        void backwards_pass(WorldSpace v1) {
            assert(samples.size() >= 1);

            size_t num_segments = samples.size() - 1;
            velocities[num_segments] = v1;
            for (size_t i = num_segments; i > 0; i--) {
                // Pass over segment i, from sample i to i-1
                PathSample sample_0 = samples[i];
                PathSample sample_1 = samples[i - 1];
                double curvature = (sample_0.curvature + sample_1.curvature) / 2.0;
                WorldSpace velocity_0 = velocities[i];

                Eigen::Matrix<double, 3, 1> direction_0_normalized = sample_0.direction;

                {
                    double direction_0_norm = direction_0_normalized.block<2, 1>(0, 0).norm();
                    if (direction_0_norm > 1e-6) {
                        direction_0_normalized /= direction_0_norm;
                    } else {
                        direction_0_normalized = sample_0.acceleration / sample_0.acceleration.block<2, 1>(0, 0).norm();
                        if (i == num_segments) {
                            direction_0_normalized = -direction_0_normalized;
                        }
                    }
                }

                // Optimize for the maximum forward acceleration, or maintain the maximum
                // velocity that we can keep without having to provide too much centripetal force.
                WorldSpace normal_vector(
                        -std::copysign(1, curvature) * direction_0_normalized(1),
                        std::copysign(1, curvature) * direction_0_normalized(0), 0);
                WorldSpace centripetal_acceleration =
                    normal_vector * velocity_0.dot(velocity_0) * curvature;

                Eigen::Matrix<double, 4, 1> wheel_acceleration;
                double ell;
                bool exact = OptimizeOnLine(
                        Eigen::Matrix<double, 4, 1>(world_to_wheel * direction_0_normalized),
                        Eigen::Matrix<double, 4, 1>(world_to_wheel * centripetal_acceleration),
                        Eigen::Matrix<double, 4, 1>(-a_max),
                        Eigen::Matrix<double, 4, 1>(a_max),
                        &wheel_acceleration,
                        &ell,
                        false);

                if (exact) {
                    WorldSpace forward_acceleration = ell * direction_0_normalized;
                    WorldSpace net_acceleration = centripetal_acceleration + forward_acceleration;

                    // Approximate the segment as a linear piece, with the forward acceleration and velocity components:
                    // dx = a*dt^2 + v*dt
                    // 0 = a*dt^2 + v*dt - dx
                    // dt = (-v +/- sqrt(v^2 + 4*a*dx))/(2a)
                    // dt must be positive, so:
                    // dt = (-v + sqrt(v^2 + 4*a*dx))/(2a)
                    double a_mag = forward_acceleration.norm();
                    double v_mag = velocity_0.norm();
                    double dx = (sample_0.position - sample_1.position).norm();
                    double dt = (-v_mag + std::sqrt(v_mag * v_mag + 4 * a_mag * dx)) / (2 * a_mag);

                    if (dt > time_deltas[i - 1]) {
                        time_deltas[i - 1] = dt;
                        velocities[i - 1] = velocity_0 + dt * net_acceleration;
                    }
                }
            }
        }

        void velocity_pass() {
            assert(samples.size() >= 1);
            size_t num_segments = samples.size() - 1;
            for (size_t i = 0; i < velocities.size(); i++) {
                // Pass over segment i, from sample i to i-1
                PathSample sample = samples[i];
                double curvature = sample.curvature;
                WorldSpace velocity_0 = velocities[i];

                Eigen::Matrix<double, 3, 1> direction_normalized = sample.direction;

                {
                    double direction_norm = direction_normalized.block<2, 1>(0, 0).norm();
                    if (direction_norm > 1e-6) {
                        direction_normalized /= direction_norm;
                    } else {
                        direction_normalized = sample.acceleration / sample.acceleration.block<2, 1>(0, 0).norm();
                        if (i == num_segments) {
                            direction_normalized = -direction_normalized;
                        }
                    }
                }

                // Optimize for the maximum forward acceleration, or maintain the maximum
                // velocity that we can keep without having to provide too much centripetal force.
                WorldSpace normal_vector(
                        -std::copysign(1, curvature) * direction_normalized(1),
                        std::copysign(1, curvature) * direction_normalized(0), 0);

                Eigen::Matrix<double, 4, 1> wheel_acceleration;
                double ell = 0;

                // Find the maximum velocity we can actually have by maximizing centripetal acceleration:
                bool exact_centripetal = OptimizeOnLine(
                        Eigen::Matrix<double, 4, 1>(world_to_wheel * normal_vector),
                        Eigen::Matrix<double, 4, 1>(Eigen::Matrix<double, 4, 1>::Zero()),
                        Eigen::Matrix<double, 4, 1>(-a_max),
                        Eigen::Matrix<double, 4, 1>(a_max),
                        &wheel_acceleration,
                        &ell,
                        false);

                assert(exact_centripetal);

                WorldSpace net_acceleration = normal_vector * ell;

                // In the case where curvature vanishes, we have unlimited speed (or v_max)
                double speed_max = v_max;
                if (std::abs(curvature) > 1e-6) {
                    // a = v^2*curvature, so v = sqrt(a / curvature)
                    speed_max = std::min(speed_max, std::sqrt(std::abs(net_acceleration.norm() / curvature)));
                }

                WorldSpace velocity = speed_max * direction_normalized;

                velocities[i] = velocity;
                if (i != 0) {
                    PathSample sample_prev = samples[i - 1];
                    double dx = (sample.position - sample_prev.position).norm();
                    double v_avg = (velocities[i - 1] + velocities[i]).norm() / 2.0;
                    time_deltas[i - 1] = dx / v_avg;
                }
            }
        }

        Eigen::Matrix<double, 4, 3> world_to_wheel;
        Eigen::Matrix<double, 4, 1> a_max;

        double v_max;

        std::vector<PathSample> samples;
        std::vector<WorldSpace> velocities;
        std::vector<double> time_deltas;
};

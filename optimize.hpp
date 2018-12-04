#include "eigen3/Eigen/Core"
#include <limits>
#include <utility>

template<int NumStates>
bool OptimizeOnLine(const Eigen::Matrix<double, NumStates, 1>& tangent,
                    const Eigen::Matrix<double, NumStates, 1>& normal,
                    const Eigen::Matrix<double, NumStates, 1>& minimum,
                    const Eigen::Matrix<double, NumStates, 1>& maximum,
                    Eigen::Matrix<double, NumStates, 1>* result,
                    double* multiplier,
                    bool do_inexact_calc = true) {
    // Maximize a subject to min <= a*tangent + normal <= max.
    // If there is no such value of a (i.e. the line formed by the tangent
    // passing through the origin does not intersect with the box generated
    // by the constraint), find the nearest corner of the constraint box.
    //
    // min - normal <= a * tangent <= max - normal
    Eigen::Matrix<double, NumStates, 1> adjusted_min = minimum - normal,
                                             adjusted_max = maximum - normal;

    double a_max = std::numeric_limits<double>::infinity();
    double a_min = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < NumStates; i++) {
        double coordinate = tangent(i);
        double a1 = adjusted_min(i) / coordinate;
        double a2 = adjusted_max(i) / coordinate;

        double a_max_candidate = std::max(a1, a2);
        double a_min_candidate = std::min(a1, a2);

        a_max = std::min(a_max_candidate, a_max);
        a_min = std::max(a_min_candidate, a_min);

        if (a_min > a_max) {
            break;
        }
    }

    if (a_min > a_max) {
        if (do_inexact_calc) {
            static_assert(NumStates <= 64);

            // Find the closest corner to the line by taking the dot product of each corner with a vector normal to the line.
            double min_dot_abs = std::numeric_limits<double>::infinity();

            Eigen::Matrix<double, NumStates, 1> best_corner;

            // Find a vector actually normal to tangent - don't just assume `normal` is valid.
            Eigen::Matrix<double, NumStates, 1> really_normal =
                normal - normal * (normal.normalized().dot(tangent.normalized()));

            assert(std::abs(really_normal.dot(tangent)) < 1e-6);

            // The bounding box doesn't intersect the line. The optimal point is at one of the corners of the box.
            // Iterate through all of the corners. Each corner is defined by a n-tuple of booleans representing which
            // bound (high or low) is active.
            for (uint64_t active_constraint_bits = 0; active_constraint_bits < (1 << NumStates); active_constraint_bits++) {
                // Construct a point at the given corner
                Eigen::Matrix<double, NumStates, 1> corner;
                for (int i = 0; i < NumStates; i++) {
                    if (active_constraint_bits >> i & 1) {
                        corner(i, 0) = adjusted_max(i, 0);
                    } else {
                        corner(i, 0) = adjusted_min(i, 0);
                    }
                }

                double dot_abs = std::abs(really_normal.dot(corner));
                if (dot_abs < min_dot_abs) {
                    min_dot_abs = dot_abs;
                    best_corner = corner;
                }
            }

            if (result) {
                *result = best_corner + normal;
            }
        }

        return false;
    } else {
        if (multiplier) {
            *multiplier = a_max;
        }
        if (result) {
            *result = normal + tangent * a_max;
        }
        return true;
    }
}

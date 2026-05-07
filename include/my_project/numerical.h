// include/my_project/numerical.h
#pragma once
#include <Eigen/Core>
#include <vector>

#include "vehicle_data.h"

namespace acv {

inline constexpr double kEps = 1e-9;

struct OscillationMode {
    std::complex<double> eigenvalue;

    double period;                 // T
    double logarithmic_decrement;  // χ
    double decay_ratio;            // K = exp(χ)
};

struct StabilityResult {
    bool is_stable;
    std::vector<std::complex<double>> eigenvalues;
    double max_real_part;

    std::vector<OscillationMode> oscillation_modes;
};

StabilityResult AnalyzeStabilitySimple(const VehicleData& data);

StabilityResult AnalyzeStabilityFull(const VehicleData& data);

}  // namespace acv

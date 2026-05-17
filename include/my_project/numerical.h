// include/my_project/numerical.h
#pragma once
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <vector>

#include "my_project/constants.h"
#include "vehicle_data.h"

namespace acv {

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

template <class MatrixType>
StabilityResult RootsProcessing(const Eigen::EigenSolver<MatrixType>& solver) {
    using namespace constants;

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigen solver failed to converge.");
    }

    StabilityResult result;
    result.max_real_part = -1e18;
    result.is_stable = true;

    auto roots = solver.eigenvalues();

    for (int i = 0; i < roots.size(); ++i) {
        std::complex<double> lambda = roots[i];
        result.eigenvalues.push_back(lambda);
        if (lambda.real() > result.max_real_part)
            result.max_real_part = lambda.real();
        if (lambda.real() > kEps) result.is_stable = false;

        if (std::abs(lambda.imag()) > kEps) {
            OscillationMode mode;
            mode.eigenvalue = lambda;
            const double sigma = lambda.real();
            const double omega = std::abs(lambda.imag());
            mode.period = 2.0 * M_PI / omega;
            mode.logarithmic_decrement = -sigma * mode.period;
            mode.decay_ratio = std::exp(mode.logarithmic_decrement);
            result.oscillation_modes.push_back(mode);
        }
    }
    return result;
}

StabilityResult AnalyzeStabilitySimple(const VehicleData& data);

StabilityResult AnalyzeStabilityFull(const VehicleData& data);

}  // namespace acv

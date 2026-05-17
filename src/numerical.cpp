// src/numerical.cpp
#include <Eigen/Core>
#include <cmath>
#include <unsupported/Eigen/Polynomials>

#include "my_project/constants.h"
#include "my_project/numerical.h"
#include "my_project/vehicle_data.h"

namespace acv {

StabilityResult AnalyzeStabilitySimple(const VehicleData& data) {
    using namespace constants;

    double compressibility = data.W0 / (kN * kPAtm);

    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

    const double Seff = 0.5 * data.L * std::pow(data.l, 2);

    A(0, 0) = (data.dQin_dp - data.dQout_dp) / compressibility;
    A(0, 1) = -data.dQout_dphi / compressibility;
    A(0, 2) = -Seff / compressibility;

    A(1, 2) = 1.0;

    const double dM_dp = 0.5 * data.L * std::pow(data.l, 2);

    A(2, 0) = dM_dp / data.I_phi;
    A(2, 1) = data.dM_dphi / data.I_phi;
    A(2, 2) = data.dM_dphidot / data.I_phi;

    Eigen::EigenSolver<Eigen::Matrix3d> solver(A);

    return RootsProcessing(solver);
}

StabilityResult AnalyzeStabilityFull(const VehicleData& data) {
    using namespace constants;

    double compressibility = data.W0 / (kN * kPAtm);

    Eigen::Matrix<double, 5, 5> A = Eigen::Matrix<double, 5, 5>::Zero();

    double Seff = 0.5 * data.L * std::pow(data.l, 2);

    A(0, 0) = (data.dQin_dp - data.dQout_dp) / compressibility;
    A(0, 1) = -data.dQout_dH / compressibility;
    A(0, 2) = -data.S / compressibility;
    A(0, 3) = -data.dQout_dphi / compressibility;
    A(0, 4) = -Seff / compressibility;

    A(1, 2) = 1.0;

    A(2, 0) = data.S / data.m;
    A(2, 2) = data.dY_dHdot / data.m;

    A(3, 4) = 1.0;

    A(4, 0) = Seff / data.I_phi;
    A(4, 3) = data.dM_dphi / data.I_phi;
    A(4, 4) = data.dM_dphidot / data.I_phi;

    // 3. Расчет собственных значений
    Eigen::EigenSolver<Eigen::Matrix<double, 5, 5>> solver(A);

    return RootsProcessing(solver);
}

}  // namespace acv

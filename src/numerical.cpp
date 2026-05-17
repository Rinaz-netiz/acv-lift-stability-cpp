// src/numerical.cpp
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <unsupported/Eigen/Polynomials>

#include "my_project/constants.h"
#include "my_project/numerical.h"
#include "my_project/vehicle_data.h"

namespace acv {

StabilityResult AnalyzeStabilitySimple(const VehicleData& data) {
    using namespace constants;

    double compressibility = data.geometry.W0 / (kN * kPAtm);
    double Seff = 0.5 * data.geometry.L * std::pow(data.geometry.l, 2);

    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

    // Уравнение давления (используем data.der)
    A(0, 0) = (data.der.dQin_dp - data.der.dQout_dp) / compressibility;
    A(0, 1) = -data.der.dQout_dphi / compressibility;
    A(0, 2) = -Seff / compressibility;

    // Кинематика
    A(1, 2) = 1.0;

    // Динамика ограждения
    A(2, 0) = Seff / data.geometry.I_phi;
    A(2, 1) = data.der.dM_dphi / data.geometry.I_phi;
    A(2, 2) = data.der.dM_dphidot / data.geometry.I_phi;

    // Решение задачи на СЗ
    Eigen::EigenSolver<Eigen::Matrix3d> solver(A);

    return RootsProcessing(solver);
}

StabilityResult AnalyzeStabilityFull(const VehicleData& data) {
    using namespace constants;

    double compressibility = data.geometry.W0 / (kN * kPAtm);
    double Seff = 0.5 * data.geometry.L * std::pow(data.geometry.l, 2);

    Eigen::Matrix<double, 5, 5> A = Eigen::Matrix<double, 5, 5>::Zero();

    // Уравнение давления
    A(0, 0) = (data.der.dQin_dp - data.der.dQout_dp) / compressibility;
    A(0, 1) = -data.der.dQout_dH / compressibility;
    A(0, 2) = -data.geometry.S / compressibility;
    A(0, 3) = -data.der.dQout_dphi / compressibility;
    A(0, 4) = -Seff / compressibility;

    // Кинематика высоты
    A(1, 2) = 1.0;

    // Динамика высоты
    A(2, 0) = data.geometry.S / data.geometry.m;
    A(2, 2) = data.der.dY_dHdot / data.geometry.m;

    // Кинематика угла
    A(3, 4) = 1.0;

    // Динамика угла
    A(4, 0) = Seff / data.geometry.I_phi;
    A(4, 3) = data.der.dM_dphi / data.geometry.I_phi;
    A(4, 4) = data.der.dM_dphidot / data.geometry.I_phi;

    Eigen::EigenSolver<Eigen::Matrix<double, 5, 5>> solver(A);

    return RootsProcessing(solver);
}

}  // namespace acv

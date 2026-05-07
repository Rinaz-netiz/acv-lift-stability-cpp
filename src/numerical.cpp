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

    // Коэффициент сжимаемости воздуха в подушке
    double compressibility = data.W0 / (kN * kPAtm);

    // 2. Формируем матрицу Якоби (Matrix of the state-space model)
    // Состояние x = [delta_p, delta_phi, delta_phi_dot]
    // Линеаризованная система: dx/dt = A * x
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

    // --- Уравнение 1: d(delta_p)/dt ---
    // (W0/n*pa) * p_dot = (dQin/dp - dQout/dp)*dp - dQout/dphi*dphi - (1/2 * L
    // * l^2)*phi_dot
    const double Seff =
        0.5 * data.L *
        std::pow(data.l, 2);  // Геометрический вклад в изменение объема

    A(0, 0) = (data.dQin_dp - data.dQout_dp) / compressibility;
    A(0, 1) = -data.dQout_dphi / compressibility;
    A(0, 2) = -Seff / compressibility;

    // --- Уравнение 2: d(delta_phi)/dt = delta_phi_dot ---
    A(1, 2) = 1.0;

    // --- Уравнение 3: d(delta_phi_dot)/dt (Уравнение моментов) ---
    // I_phi * phi_ddot = (L * l^2 / 2) * dp + dM/dphi * dphi + dM/dphidot *
    // dphi_dot
    const double dM_dp =
        0.5 * data.L * std::pow(data.l, 2);  // Момент от давления

    A(2, 0) = dM_dp / data.I_phi;
    A(2, 1) = data.dM_dphi / data.I_phi;
    A(2, 2) = data.dM_dphidot / data.I_phi;

    // 3. Решаем задачу на собственные значения
    Eigen::EigenSolver<Eigen::Matrix3d> solver(A);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigen solver failed to converge.");
    }

    StabilityResult result;
    result.max_real_part = -1e18;
    result.is_stable = true;

    auto evals = solver.eigenvalues();
    for (int i = 0; i < evals.size(); ++i) {
        std::complex<double> lambda = evals[i];
        result.eigenvalues.push_back(lambda);

        if (lambda.real() > result.max_real_part) {
            result.max_real_part = lambda.real();
        }

        // Условие устойчивости: Вещественная часть < 0
        if (lambda.real() > kEps) {
            result.is_stable = false;
        }

        if (lambda.imag() > kEps) {
            OscillationMode mode;

            mode.eigenvalue = lambda;

            const double sigma = lambda.real();
            const double omega = std::abs(lambda.imag());

            // T = 2π / ω
            mode.period = 2.0 * M_PI / omega;

            // χ = -σT
            mode.logarithmic_decrement = -sigma * mode.period;

            // K = exp(χ)
            mode.decay_ratio = std::exp(mode.logarithmic_decrement);

            result.oscillation_modes.push_back(mode);
        }
    }

    return result;
}

StabilityResult AnalyzeStabilityFull(const VehicleData& data) {
    using namespace constants;

    // Коэффициент сжимаемости объема воздуха в подушке
    double compressibility = data.W0 / (kN * kPAtm);

    // 2. Инициализация матрицы Якоби 5x5
    // Состояния: 0:dp, 1:dH, 2:dH_dot, 3:dphi, 4:dphi_dot
    Eigen::Matrix<double, 5, 5> A = Eigen::Matrix<double, 5, 5>::Zero();

    // Геометрический вклад гибкого ограждения в изменение объема (S_eff)
    // Согласно статье, это L * l^2 / 2
    double Seff = 0.5 * data.L * std::pow(data.l, 2);

    // --- Уравнение 1: d(delta_p)/dt (Уравнение неразрывности) ---
    // Формула: compressibility * p_dot = (dQin/dp - dQout/dp)*dp - dQout/dH*dH
    // - S*H_dot - dQout/dphi*dphi - Seff*phi_dot
    A(0, 0) = (data.dQin_dp - data.dQout_dp) / compressibility;
    A(0, 1) = -data.dQout_dH / compressibility;
    A(0, 2) = -data.S / compressibility;
    A(0, 3) = -data.dQout_dphi / compressibility;
    A(0, 4) = -Seff / compressibility;

    // --- Уравнение 2: d(delta_H)/dt = delta_H_dot ---
    A(1, 2) = 1.0;

    // --- Уравнение 3: d(delta_H_dot)/dt (Вертикальная динамика корпуса) ---
    // m * H_ddot = S * dp + dY/dHdot * H_dot
    A(2, 0) = data.S / data.m;
    A(2, 2) = data.dY_dHdot / data.m;

    // --- Уравнение 4: d(delta_phi)/dt = delta_phi_dot ---
    A(3, 4) = 1.0;

    // --- Уравнение 5: d(delta_phi_dot)/dt (Динамика ограждения) ---
    // I_phi * phi_ddot = Seff * dp + dM/dphi * dphi + dM/dphidot * phi_dot
    A(4, 0) = Seff / data.I_phi;
    A(4, 3) = data.dM_dphi / data.I_phi;
    A(4, 4) = data.dM_dphidot / data.I_phi;

    // 3. Расчет собственных значений
    Eigen::EigenSolver<Eigen::Matrix<double, 5, 5>> solver(A);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigen solver failed to converge.");
    }

    StabilityResult result;
    result.max_real_part = -1e18;
    result.is_stable = true;

    auto evals = solver.eigenvalues();
    for (int i = 0; i < evals.size(); ++i) {
        std::complex<double> lambda = evals[i];
        result.eigenvalues.push_back(lambda);

        if (lambda.real() > result.max_real_part) {
            result.max_real_part = lambda.real();
        }

        // Система устойчива, если вещественные части всех СЗ строго меньше нуля
        if (lambda.real() > kEps) {
            result.is_stable = false;
        }

        if (lambda.imag() > kEps) {
            OscillationMode mode;

            mode.eigenvalue = lambda;

            const double sigma = lambda.real();
            const double omega = std::abs(lambda.imag());

            // T = 2π / ω
            mode.period = 2.0 * M_PI / omega;

            // χ = -σT
            mode.logarithmic_decrement = -sigma * mode.period;

            // K = exp(χ)
            mode.decay_ratio = std::exp(mode.logarithmic_decrement);

            result.oscillation_modes.push_back(mode);
        }
    }

    return result;
}

}  // namespace acv

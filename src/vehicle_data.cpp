// src/vehicle_data.cpp
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>

#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {

// ============================================================
// BlowerCurve
// ============================================================

// Линейная интерполяция по таблице
static double LinearInterp(const std::vector<double>& xs,
                           const std::vector<double>& ys, double x) {
    if (xs.size() < 2)
        throw std::runtime_error("Interpolation table too small");

    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();

    // Найти интервал
    size_t i = 1;
    while (i < xs.size() - 1 && xs[i] < x) ++i;

    double t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
    return ys[i - 1] + t * (ys[i] - ys[i - 1]);
}

static double LinearInterpDerivative(const std::vector<double>& xs,
                                     const std::vector<double>& ys, double x) {
    if (xs.size() < 2)
        throw std::runtime_error("Interpolation table too small");

    if (x <= xs.front()) return (ys[1] - ys[0]) / (xs[1] - xs[0]);
    if (x >= xs.back()) {
        size_t n = xs.size();
        return (ys[n - 1] - ys[n - 2]) / (xs[n - 1] - xs[n - 2]);
    }

    size_t i = 1;
    while (i < xs.size() - 1 && xs[i] < x) ++i;

    return (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
}

double BlowerCurve::Flow(double p) const {
    if (!p_table.empty()) {
        return LinearInterp(p_table, Q_table, p);
    }
    // Линейная модель: Q = Q_max - k_fan * p
    return Q_max - k_fan * p;
}

double BlowerCurve::Derivative(double p) const {
    if (!p_table.empty()) {
        return LinearInterpDerivative(p_table, Q_table, p);
    }
    return -k_fan;  // dQ/dp = -k_fan < 0
}

// ============================================================
// SealMomentCurve
// M(φ) — момент сопротивления ограждения (< 0 при φ > 0)
// В таблице хранится M/L, итоговый момент = (M/L) * L
// ============================================================

double SealMomentCurve::Moment(double phi, double L) const {
    if (phi_table.empty())
        throw std::runtime_error("SealMomentCurve: table is empty");

    double M_per_L = LinearInterp(phi_table, M_per_L_table, phi);
    return M_per_L * L;
}

double SealMomentCurve::Derivative(double phi, double L) const {
    if (phi_table.empty())
        throw std::runtime_error("SealMomentCurve: table is empty");

    double dM_per_L_dphi =
        LinearInterpDerivative(phi_table, M_per_L_table, phi);
    return dM_per_L_dphi * L;  // < 0
}

// ============================================================
// Eq. (24): зазор под ограждением
// h_gap = H - H_init + l*(1 - cos(φ))
// ============================================================
static double ComputeGap(const VehicleGeometry& g, double H, double phi) {
    return (H - g.H_init) + g.l * (1.0 - std::cos(phi));
}

// Eq. (11): расход утечки
// Q_out = χ * sqrt(2p/ρ) * S_gap = χ * sqrt(2p/ρ) * L * h_gap
static double ComputeQout(double p, double h_gap, double L) {
    using namespace constants;
    if (h_gap <= 0.0) return 0.0;
    return kChi * std::sqrt(2.0 * p / kRhoAir) * L * h_gap;
}

// ============================================================
// Eq. (33): система уравнений равновесия
//
//   (1)  p * S = m * g
//   (2)  Q_in(p) = χ*sqrt(2p/ρ)*L*(H - H_init + l*(1-cos(φ)))
//   (3)  p * L * l²/2 = -M(φ, 0)
//
// Из (1): p0 = m*g/S  — явно
// Из (3): φ0 — нелинейное уравнение (M(φ) нелинейна)
// Из (2): H0 — явно при известных p0, φ0
// ============================================================

void VehicleData::ComputeEquilibrium() {
    using namespace constants;

    const double m = geometry.m;
    const double g = kG;
    const double S = geometry.S;
    const double L = geometry.L;
    const double l = geometry.l;

    // ---- Шаг 1: p0 из Eq. (33, строка 1) ----
    // p0 * S = m * g
    eq.p0 = (m * g) / S;

    if (eq.p0 <= 0.0) throw std::runtime_error("Equilibrium pressure p0 <= 0");

    // ---- Шаг 2: φ0 из Eq. (33, строка 3) ----
    // p0 * L * l²/2 + M(φ0) = 0
    // M(φ0) = -p0 * L * l²/2
    //
    // Решаем методом бисекции: f(φ) = p0*L*l²/2 + M(φ) = 0
    // M(φ) < 0, строго убывает → f монотонна

    const double M_target = -(eq.p0 * L * l * l / 2.0);
    // Ищем φ такой, что M(φ)*L = M_target

    // Проверяем диапазон таблицы
    double phi_lo = seal_moment.phi_table.front();
    double phi_hi = seal_moment.phi_table.back();

    auto f = [&](double phi) -> double {
        // M(phi) + p0*L*l²/2 = 0
        return seal_moment.Moment(phi, L) + (eq.p0 * L * l * l / 2.0);
    };

    double f_lo = f(phi_lo);
    double f_hi = f(phi_hi);

    if (f_lo * f_hi > 0.0) {
        throw std::runtime_error(
            "Cannot find equilibrium angle: M(phi) curve does not "
            "intersect required value M_target = " +
            std::to_string(M_target) +
            " N·m/m.\n"
            "  M(phi_lo=" +
            std::to_string(phi_lo) +
            ") = " + std::to_string(seal_moment.Moment(phi_lo, L)) +
            "\n"
            "  M(phi_hi=" +
            std::to_string(phi_hi) +
            ") = " + std::to_string(seal_moment.Moment(phi_hi, L)) +
            "\n"
            "  Required: " +
            std::to_string(M_target));
    }

    // Бисекция — надёжна, т.к. M(φ) монотонна
    const int bisect_max = 100;
    const double bisect_tol = 1e-10;

    for (int i = 0; i < bisect_max; ++i) {
        double phi_mid = 0.5 * (phi_lo + phi_hi);
        double f_mid = f(phi_mid);

        if (std::abs(f_mid) < bisect_tol || (phi_hi - phi_lo) < bisect_tol) {
            eq.phi0 = phi_mid;
            break;
        }

        if (f_lo * f_mid < 0.0) {
            phi_hi = phi_mid;
            f_hi = f_mid;
        } else {
            phi_lo = phi_mid;
            f_lo = f_mid;
        }

        if (i == bisect_max - 1) eq.phi0 = 0.5 * (phi_lo + phi_hi);
    }

    // ---- Шаг 3: H0 из Eq. (33, строка 2) ----
    // Q_in(p0) = χ*sqrt(2p0/ρ)*L*(H0 - H_init + l*(1 - cos(φ0)))
    // => h_gap0 = Q_in(p0) / (χ*sqrt(2p0/ρ)*L)
    // => H0 = H_init + h_gap0 - l*(1 - cos(φ0))

    double Q_in0 = blower.Flow(eq.p0);

    if (Q_in0 <= 0.0) {
        throw std::runtime_error(
            "Blower flow Q_in(p0) <= 0 at p0 = " + std::to_string(eq.p0) +
            " Pa. "
            "Increase Q_max or decrease k_fan.");
    }

    double velocity = kChi * std::sqrt(2.0 * eq.p0 / kRhoAir);

    eq.h_gap0 = Q_in0 / (velocity * geometry.L);

    if (eq.h_gap0 <= 0.0) {
        throw std::runtime_error("Equilibrium gap h_gap0 <= 0");
    }

    eq.H0 =
        geometry.H_init + eq.h_gap0 - geometry.l * (1.0 - std::cos(eq.phi0));

    eq.Q0 = Q_in0;  // В равновесии Q_in = Q_out = Q0

    // ---- Проверка баланса ----
    {
        double h_check = ComputeGap(geometry, eq.H0, eq.phi0);
        double Q_out_check = ComputeQout(eq.p0, h_check, geometry.L);
        double err = std::abs(Q_out_check - eq.Q0) / eq.Q0;
        if (err > 1e-6) {
            throw std::runtime_error(
                "Flow balance check failed: relative error = " +
                std::to_string(err));
        }
    }
}

// ============================================================
// Производные в точке равновесия — Eq. (25)
// ============================================================

void VehicleData::ComputeDerivatives() {
    using namespace constants;

    // dQout/dp|0 = Q0 / (2*p0)   [Eq. 25]
    der.dQout_dp = eq.Q0 / (2.0 * eq.p0);

    // dQout/dH|0 = χ * sqrt(2p0/ρ) * L   [Eq. 25]
    der.dQout_dH = kChi * std::sqrt(2.0 * eq.p0 / kRhoAir) * geometry.L;

    // dQout/dφ|0 = χ * sqrt(2p0/ρ) * L * l * sin(φ0)   [Eq. 25]
    der.dQout_dphi = kChi * std::sqrt(2.0 * eq.p0 / kRhoAir) * geometry.L *
                     geometry.l * std::sin(eq.phi0);

    // dQin/dp|0 — из характеристики нагнетателя [Eq. 12]
    der.dQin_dp = blower.Derivative(eq.p0);  // < 0

    // dM/dφ|0 — из кривой момента [Section 2.2.1, Fig. 6]
    der.dM_dphi = seal_moment.Derivative(eq.phi0, geometry.L);  // < 0

    // dM/dφ̇|0 — входное данное [Section 2.2.2, Table 1]
    der.dM_dphidot = dM_dphidot_input;  // < 0

    // dY/dḢ|0 = D(h_gap0) * ρ * (-Q_in/S) * S = -D * ρ * Q0
    // Eq. (13): p_damp = D(S_gap)*ρ*(-Q_in/S * H_dot + ...)*S
    // Линеаризация: dY/dH_dot = -D * ρ * Q0
    der.dY_dHdot = -D_damping * kRhoAir * eq.Q0;  // < 0
}

// ============================================================
// Валидация
// ============================================================

void VehicleData::Validate() const {
    // Геометрия
    if (geometry.m <= 0) throw std::runtime_error("m <= 0");
    if (geometry.L <= 0) throw std::runtime_error("L <= 0");
    if (geometry.l <= 0) throw std::runtime_error("l <= 0");
    if (geometry.S <= 0) throw std::runtime_error("S <= 0");
    if (geometry.W0 <= 0) throw std::runtime_error("W0 <= 0");
    if (geometry.I_phi <= 0) throw std::runtime_error("I_phi <= 0");

    // Нагнетатель
    if (blower.p_table.empty()) {
        if (blower.Q_max <= 0) throw std::runtime_error("Q_max <= 0");
        if (blower.k_fan <= 0) throw std::runtime_error("k_fan <= 0");
    } else {
        if (blower.p_table.size() != blower.Q_table.size())
            throw std::runtime_error("Blower table sizes mismatch");
        if (blower.p_table.size() < 2)
            throw std::runtime_error("Blower table too small");
    }

    // Таблица момента
    if (seal_moment.phi_table.empty())
        throw std::runtime_error("Seal moment table is empty");
    if (seal_moment.phi_table.size() != seal_moment.M_per_L_table.size())
        throw std::runtime_error("Seal moment table sizes mismatch");

    // Демпфирование
    if (dM_dphidot_input >= 0.0)
        throw std::runtime_error("dM_dphidot must be < 0");
    if (D_damping < 0.0) throw std::runtime_error("D_damping must be >= 0");

    // Если равновесие вычислено
    if (eq.p0 > 0.0) {
        if (eq.Q0 <= 0.0) throw std::runtime_error("Q0 <= 0");
        if (eq.h_gap0 <= 0.0) throw std::runtime_error("h_gap0 <= 0");
        if (eq.phi0 < 0.0) throw std::runtime_error("phi0 < 0");

        // Eq. (3) — необходимое условие устойчивости подушки
        double crit = der.dQin_dp - der.dQout_dp;
        if (crit >= 0.0)
            throw std::runtime_error(
                "Necessary stability criterion violated: "
                "dQin/dp - dQout/dp = " +
                std::to_string(crit) + " >= 0");

        if (der.dM_dphi >= 0.0) throw std::runtime_error("dM/dphi must be < 0");
        if (der.dY_dHdot >= 0.0)
            throw std::runtime_error("dY/dHdot must be < 0");
    }
}

void VehicleData::Init() {
    Validate();            // Проверка входных данных
    ComputeEquilibrium();  // Eq. (33) → p0, φ0, H0, Q0
    ComputeDerivatives();  // Eq. (25) → производные
    Validate();            // Финальная проверка
}

// ============================================================
// JSON
// ============================================================

VehicleData LoadVehicleFromJson(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    nlohmann::json j;
    file >> j;

    VehicleData v;
    v.name = j["name"];

    // Геометрия
    auto& g = j["geometry"];
    v.geometry.m = g["m"];
    v.geometry.L = g["L"];
    v.geometry.l = g["l"];
    v.geometry.S = g["S"];
    v.geometry.W0 = g["W0"];
    v.geometry.I_phi = g["I_phi"];
    v.geometry.H_init = g["H_init"];

    // Нагнетатель
    auto& bl = j["blower"];
    if (bl.contains("p_table")) {
        // Табличная характеристика
        v.blower.p_table = bl["p_table"].get<std::vector<double>>();
        v.blower.Q_table = bl["Q_table"].get<std::vector<double>>();
    } else {
        // Линейная модель
        v.blower.Q_max = bl["Q_max"];
        v.blower.k_fan = bl["k_fan"];
    }

    // Момент ограждения — таблица (φ, M/L)
    auto& sm = j["seal_moment"];
    v.seal_moment.phi_table = sm["phi"].get<std::vector<double>>();
    v.seal_moment.M_per_L_table = sm["M_per_L"].get<std::vector<double>>();

    // Демпфирование
    v.dM_dphidot_input = j["dM_dphidot"];  // < 0, из CFD / Table 1
    v.D_damping = j["D_damping"];          // > 0, из Fig. 4

    v.Init();
    return v;
}

void SaveVehicleToJson(const VehicleData& v, const std::string& filename) {
    nlohmann::json j;
    j["name"] = v.name;

    j["geometry"]["m"] = v.geometry.m;
    j["geometry"]["L"] = v.geometry.L;
    j["geometry"]["l"] = v.geometry.l;
    j["geometry"]["S"] = v.geometry.S;
    j["geometry"]["W0"] = v.geometry.W0;
    j["geometry"]["I_phi"] = v.geometry.I_phi;
    j["geometry"]["H_init"] = v.geometry.H_init;

    if (!v.blower.p_table.empty()) {
        j["blower"]["p_table"] = v.blower.p_table;
        j["blower"]["Q_table"] = v.blower.Q_table;
    } else {
        j["blower"]["Q_max"] = v.blower.Q_max;
        j["blower"]["k_fan"] = v.blower.k_fan;
    }

    j["seal_moment"]["phi"] = v.seal_moment.phi_table;
    j["seal_moment"]["M_per_L"] = v.seal_moment.M_per_L_table;

    j["dM_dphidot"] = v.dM_dphidot_input;
    j["D_damping"] = v.D_damping;

    // Вычисленные значения — для информации
    j["equilibrium"]["p0"] = v.eq.p0;
    j["equilibrium"]["phi0"] = v.eq.phi0;
    j["equilibrium"]["H0"] = v.eq.H0;
    j["equilibrium"]["Q0"] = v.eq.Q0;
    j["equilibrium"]["h_gap0"] = v.eq.h_gap0;

    j["derivatives"]["dQin_dp"] = v.der.dQin_dp;
    j["derivatives"]["dQout_dp"] = v.der.dQout_dp;
    j["derivatives"]["dQout_dH"] = v.der.dQout_dH;
    j["derivatives"]["dQout_dphi"] = v.der.dQout_dphi;
    j["derivatives"]["dM_dphi"] = v.der.dM_dphi;
    j["derivatives"]["dM_dphidot"] = v.der.dM_dphidot;
    j["derivatives"]["dY_dHdot"] = v.der.dY_dHdot;

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error("Cannot write file: " + filename);
    out << j.dump(4);
}

}  // namespace acv

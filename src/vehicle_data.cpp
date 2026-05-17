// src/vehicle_data.cpp

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iterator>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <utility>
#include <vector>

#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {

// Linear interpolation
static size_t FindInterval(const std::vector<double>& xs, double x) {
    auto it = std::lower_bound(xs.begin(), xs.end(), x);

    if (it == xs.begin()) {
        return 1;
    }
    if (it == xs.end()) {
        return xs.size() - 1;
    }

    return std::distance(xs.begin(), it);
}

static double LinearInterp(const std::vector<double>& xs,
                           const std::vector<double>& ys, double x) {
    const size_t n = xs.size();
    if (n < 2) {
        throw std::runtime_error("Interpolation table too small");
    }

    if (x <= xs.front()) {
        return ys.front();
    }
    if (x >= xs.back()) {
        return ys.back();
    }

    size_t i = FindInterval(xs, x);

    double dx = xs[i] - xs[i - 1];
    if (dx == 0.0) {
        return ys[i - 1];
    }

    double t = (x - xs[i - 1]) / dx;
    return ys[i - 1] + t * (ys[i] - ys[i - 1]);
}

static double LinearInterpDerivative(const std::vector<double>& xs,
                                     const std::vector<double>& ys, double x) {
    const size_t n = xs.size();
    if (n < 2) throw std::runtime_error("Interpolation table too small");

    size_t i;
    if (x <= xs.front()) {
        i = 1;
    } else if (x >= xs.back()) {
        i = n - 1;
    } else {
        i = FindInterval(xs, x);
    }

    double dx = xs[i] - xs[i - 1];
    if (dx == 0.0) {
        return 0.0;
    }

    return (ys[i] - ys[i - 1]) / dx;
}

// BlowerCurve
double BlowerCurve::Flow(double p) const {
    if (!p_table.empty()) {
        return LinearInterp(p_table, Q_table, p);
    }
    return Q_max - k_fan * p;
}

double BlowerCurve::Derivative(double p) const {
    if (!p_table.empty()) {
        return LinearInterpDerivative(p_table, Q_table, p);
    }
    return -k_fan;
}

double SealDampingCurve::Value(double h) const {
    return LinearInterp(h_gap_table, dM_dphidot_table, h);
}

double CushionDampingCurve::Value(double r) const {
    return LinearInterp(Sgap_over_S_table, D_table, r);
}

// SealMomentCurve
double SealMomentCurve::Moment(double phi, double L) const {
    if (phi_table.empty()) {
        throw std::runtime_error("SealMomentCurve: table is empty");
    }

    double M_per_L = LinearInterp(phi_table, M_per_L_table, phi);
    return M_per_L * L;
}

double SealMomentCurve::Derivative(double phi, double L) const {
    if (phi_table.empty())
        throw std::runtime_error("SealMomentCurve: table is empty");

    double dM_per_L_dphi =
        LinearInterpDerivative(phi_table, M_per_L_table, phi);
    return dM_per_L_dphi * L;
}

// Eq. (24): gap under the seal
// h_gap = (H - H_init) + l*(1 - cos(φ))
static double ComputeGap(const VehicleGeometry& g, double H, double phi) {
    return (H - g.H_init) + g.l * (1.0 - std::cos(phi));
}

// Eq. (11): outflow rate
// Q_out = χ * sqrt(2p/ρ) * L * h_gap
static double ComputeQout(double p, double h_gap, double L) {
    using namespace constants;
    if (h_gap <= 0.0) return 0.0;
    return kChi * std::sqrt(2.0 * p / kRhoAir) * L * h_gap;
}

// Eq. (33): equilibrium conditions
//
//   (1)  p * S = m * g
//   (2)  Q_in(p) = χ*sqrt(2p/ρ)*L*(H - H_init + l*(1-cos(φ)))
//   (3)  p * L * l²/2 = -M(φ, 0)
void VehicleData::ComputeEquilibrium() {
    using namespace constants;

    const double m = geometry.m;
    const double g = kG;
    const double S = geometry.S;
    const double L = geometry.L;
    const double l = geometry.l;

    // ---- Step 1: p0 from Eq. (33, line 1) ----
    eq.p0 = (m * g) / S;

    if (eq.p0 <= 0.0) {
        throw std::runtime_error("Equilibrium pressure p0 <= 0");
    }

    // ---- Step 2: φ0 from Eq. (33, line 3) ----
    // p0 * L * l²/2 + M(φ0, 0) = 0
    // => M(φ0, 0) = -p0 * L * l²/2

    const double M_target = -(eq.p0 * L * l * l / 2.0);

    double phi_lo = seal_moment.phi_table.front();
    double phi_hi = seal_moment.phi_table.back();

    auto f = [&](double phi) -> double {
        return seal_moment.Moment(phi, L) - M_target;
    };

    double f_lo = f(phi_lo);
    double f_hi = f(phi_hi);

    if (f_lo * f_hi > 0.0) {
        throw std::runtime_error(
            "Cannot find equilibrium angle: M(phi) curve does not "
            "intersect required value M_target = " +
            std::to_string(M_target) +
            " N·m.\n"
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

    const int bisect_max = 200;
    const double bisect_tol = 1e-12;

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

    // ---- Step 3: H0 from Eq. (33, line 2) ----
    // Q_in(p0) = χ*sqrt(2p0/ρ)*L*(H0 - H_init + l*(1 - cos(φ0)))
    // => h_gap0 = Q_in(p0) / (χ*sqrt(2p0/ρ)*L)
    // => H0 = H_init + h_gap0 - l*(1 - cos(φ0))

    double Q_in0 = blower.Flow(eq.p0);

    if (Q_in0 <= 0.0) {
        throw std::runtime_error("Blower flow Q_in(p0) <= 0 at p0 = " +
                                 std::to_string(eq.p0) + " Pa.");
    }

    // χ * sqrt(2*p0/ρ)
    double vel = kChi * std::sqrt(2.0 * eq.p0 / kRhoAir);

    eq.h_gap0 = Q_in0 / (vel * geometry.L);

    if (eq.h_gap0 <= 0.0)
        throw std::runtime_error("Equilibrium gap h_gap0 <= 0");

    eq.H0 =
        geometry.H_init + eq.h_gap0 - geometry.l * (1.0 - std::cos(eq.phi0));

    eq.Q0 = Q_in0;  // At equilibrium: Q_in = Q_out = Q0

    // ---- Flow balance check ----

    double h_check = ComputeGap(geometry, eq.H0, eq.phi0);
    double Q_out_check = ComputeQout(eq.p0, h_check, geometry.L);
    double err = std::abs(Q_out_check - eq.Q0) / eq.Q0;
    if (err > 1e-6) {
        throw std::runtime_error(
            "Flow balance check failed: relative error = " +
            std::to_string(err));
    }
}

// Derivatives at equilibrium — Eq. (25) and related
void VehicleData::ComputeDerivatives() {
    using namespace constants;

    // ----------------------------------------------------------
    // Eq. (25): derivatives of Q_out
    //
    //   ∂Qout/∂p  |₀ = Q₀ / (2·p₀)
    //   ∂Qout/∂H  |₀ = χ₀·√(2p₀/ρ)·L
    //   ∂Qout/∂φ  |₀ = χ₀·√(2p₀/ρ)·L·l·sin(φ₀)
    // ----------------------------------------------------------

    const double chi_vel = kChi * std::sqrt(2.0 * eq.p0 / kRhoAir);
    // χ₀·√(2p₀/ρ) — used repeatedly below

    der.dQout_dp = eq.Q0 / (2.0 * eq.p0);

    der.dQout_dH = chi_vel * geometry.L;

    der.dQout_dphi = chi_vel * geometry.L * geometry.l * std::sin(eq.phi0);

    // ∂Qin/∂p|₀ — from blower curve
    der.dQin_dp = blower.Derivative(eq.p0);  // < 0

    // ∂M/∂φ|₀ — from seal moment curve
    der.dM_dphi = seal_moment.Derivative(eq.phi0, geometry.L);  // < 0

    // ----------------------------------------------------------
    // ∂M/∂φ̇|₀ — from seal damping curve (Fig. 9) or constant
    //
    // Paper Fig. 9: the curve is given as M(h_gap) in N·m
    // for a unit angular velocity, so it is already ∂M/∂φ̇.
    // ----------------------------------------------------------
    if (!seal_damping_curve.h_gap_table.empty()) {
        der.dM_dphidot = seal_damping_curve.Value(eq.h_gap0);
    } else {
        der.dM_dphidot = dM_dphidot_input;  // < 0
    }

    // ----------------------------------------------------------
    // ∂Y/∂Ḣ|₀  — linearisation of Eq. (13)
    //
    // p_damp = D(S_gap)·ρ·(-Q_in/S·Ḣ + ½·Ḣ²)·S
    // Y      = p_damp · S   (damping force on hull)
    //
    // Linearising around Ḣ = 0:
    //   ∂Y/∂Ḣ|₀ = D(S_gap₀)·ρ·(-Q₀/S)·S
    //            = -D₀·ρ·Q₀          [< 0]
    //
    // D₀ is read from Fig. 4 at S_gap₀/S.
    // ----------------------------------------------------------
    double D0 = D_damping;  // fallback constant

    if (!cushion_damping_curve.Sgap_over_S_table.empty()) {
        const double S_gap0 =
            geometry.L * eq.h_gap0;                // S_gap = L·h_gap  [Eq.24]
        const double ratio = S_gap0 / geometry.S;  // S_gap/S  [Fig. 4 x-axis]
        D0 = cushion_damping_curve.Value(ratio);
    }

    der.dY_dHdot = -D0 * kRhoAir * eq.Q0;  // < 0
}

void SortingInputedCurve(std::vector<double>& arguments,
                         std::vector<double>& values) {
    if (values.size() != arguments.size()) {
        throw std::runtime_error("table sizes mismatch");
    }

    std::vector<std::pair<double, double>> vec(arguments.size());

    for (size_t ind = 0; ind < vec.size(); ++ind) {
        vec[ind].first = arguments[ind];
        vec[ind].second = values[ind];
    }

    std::sort(vec.begin(), vec.end());

    for (size_t ind = 0; ind < vec.size(); ++ind) {
        arguments[ind] = vec[ind].first;
        values[ind] = vec[ind].second;
    }
}

// Validation
void VehicleData::Validate() const {
    // Geometry
    if (geometry.m <= 0) {
        throw std::runtime_error("m <= 0");
    }
    if (geometry.L <= 0) {
        throw std::runtime_error("L <= 0");
    }
    if (geometry.l <= 0) {
        throw std::runtime_error("l <= 0");
    }
    if (geometry.S <= 0) {
        throw std::runtime_error("S <= 0");
    }
    if (geometry.W0 <= 0) {
        throw std::runtime_error("W0 <= 0");
    }
    if (geometry.I_phi <= 0) {
        throw std::runtime_error("I_phi <= 0");
    }

    // Blower
    if (blower.p_table.empty()) {
        if (blower.Q_max <= 0) throw std::runtime_error("Q_max <= 0");
        if (blower.k_fan <= 0) throw std::runtime_error("k_fan <= 0");
    } else {
        if (blower.p_table.size() != blower.Q_table.size())
            throw std::runtime_error("Blower table sizes mismatch");
        if (blower.p_table.size() < 2)
            throw std::runtime_error("Blower table too small");
    }

    // Seal moment table
    if (seal_moment.phi_table.empty()) {
        throw std::runtime_error("Seal moment table is empty");
    }
    if (seal_moment.phi_table.size() != seal_moment.M_per_L_table.size()) {
        throw std::runtime_error("Seal moment table sizes mismatch");
    }

    // Seal damping
    const bool has_seal_curve = !seal_damping_curve.h_gap_table.empty() ||
                                !seal_damping_curve.dM_dphidot_table.empty();

    if (has_seal_curve) {
        if (seal_damping_curve.h_gap_table.size() !=
            seal_damping_curve.dM_dphidot_table.size())
            throw std::runtime_error("Seal damping curve sizes mismatch");
        if (seal_damping_curve.h_gap_table.size() < 2)
            throw std::runtime_error("Seal damping curve too small");
    } else {
        if (dM_dphidot_input >= 0.0)
            throw std::runtime_error("dM_dphidot must be < 0");
    }

    // Cushion damping
    const bool has_cushion_curve =
        !cushion_damping_curve.Sgap_over_S_table.empty() ||
        !cushion_damping_curve.D_table.empty();

    if (has_cushion_curve) {
        if (cushion_damping_curve.Sgap_over_S_table.size() !=
            cushion_damping_curve.D_table.size()) {
            throw std::runtime_error("Cushion damping curve sizes mismatch");
        }
        if (cushion_damping_curve.Sgap_over_S_table.size() < 2) {
            throw std::runtime_error("Cushion damping curve too small");
        }
    } else {
        if (D_damping < 0.0) {
            throw std::runtime_error("D_damping must be >= 0");
        }
    }

    // Post-equilibrium checks
    if (eq.p0 > 0.0) {
        if (eq.Q0 <= 0.0) {
            throw std::runtime_error("Q0 <= 0");
        }
        if (eq.h_gap0 <= 0.0) {
            throw std::runtime_error("h_gap0 <= 0");
        }
        if (eq.phi0 < 0.0) {
            throw std::runtime_error("phi0 < 0");
        }

        // Necessary stability condition: dQin/dp - dQout/dp < 0
        double crit = der.dQin_dp - der.dQout_dp;
        if (crit >= 0.0)
            throw std::runtime_error(
                "Necessary stability criterion violated: "
                "dQin/dp - dQout/dp = " +
                std::to_string(crit) + " >= 0");

        if (der.dM_dphi >= 0.0) {
            throw std::runtime_error("dM/dphi must be < 0");
        }
        if (der.dY_dHdot >= 0.0) {
            throw std::runtime_error("dY/dHdot must be < 0");
        }
    }
}

void VehicleData::Init() {
    Validate();
    ComputeEquilibrium();
    ComputeDerivatives();
    Validate();
}

// JSON I/O
VehicleData LoadVehicleFromJson(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open file: " + filename);

    nlohmann::json j;
    file >> j;

    VehicleData v;
    v.name = j["name"];

    // Geometry
    auto& g = j["geometry"];
    v.geometry.m = g["m"];
    v.geometry.L = g["L"];
    v.geometry.l = g["l"];
    v.geometry.S = g["S"];
    v.geometry.W0 = g["W0"];
    v.geometry.I_phi = g["I_phi"];
    v.geometry.H_init = g["H_init"];

    // Blower
    auto& bl = j["blower"];
    if (bl.contains("p_table")) {
        v.blower.p_table = bl["p_table"].get<std::vector<double>>();
        v.blower.Q_table = bl["Q_table"].get<std::vector<double>>();
        SortingInputedCurve(v.blower.p_table, v.blower.Q_table);
    } else {
        v.blower.Q_max = bl["Q_max"];
        v.blower.k_fan = bl["k_fan"];
    }

    // Seal moment table (φ, M/L)
    auto& sm = j["seal_moment"];
    v.seal_moment.phi_table = sm["phi"].get<std::vector<double>>();
    v.seal_moment.M_per_L_table = sm["M_per_L"].get<std::vector<double>>();

    // Seal damping curve  [Fig. 9: h_gap → ∂M/∂φ̇]
    if (j.contains("seal_damping_curve")) {
        auto& sd = j["seal_damping_curve"];
        v.seal_damping_curve.h_gap_table =
            sd["h_gap_table"].get<std::vector<double>>();
        v.seal_damping_curve.dM_dphidot_table =
            sd["dM_dphidot_table"].get<std::vector<double>>();
        SortingInputedCurve(v.seal_damping_curve.h_gap_table,
                            v.seal_damping_curve.dM_dphidot_table);
    }
    if (j.contains("dM_dphidot")) {
        v.dM_dphidot_input = j["dM_dphidot"];
    }

    // Cushion damping curve  [Fig. 4: S_gap/S → D]
    if (j.contains("cushion_damping_curve")) {
        auto& cd = j["cushion_damping_curve"];
        v.cushion_damping_curve.Sgap_over_S_table =
            cd["Sgap_over_S_table"].get<std::vector<double>>();
        v.cushion_damping_curve.D_table =
            cd["D_table"].get<std::vector<double>>();
        SortingInputedCurve(v.cushion_damping_curve.Sgap_over_S_table,
                            v.cushion_damping_curve.D_table);
    }
    if (j.contains("D_damping")) {
        v.D_damping = j["D_damping"];
    }

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

    if (!v.seal_damping_curve.h_gap_table.empty()) {
        j["seal_damping_curve"]["h_gap_table"] =
            v.seal_damping_curve.h_gap_table;
        j["seal_damping_curve"]["dM_dphidot_table"] =
            v.seal_damping_curve.dM_dphidot_table;
    } else {
        j["dM_dphidot"] = v.dM_dphidot_input;
    }

    if (!v.cushion_damping_curve.Sgap_over_S_table.empty()) {
        j["cushion_damping_curve"]["Sgap_over_S_table"] =
            v.cushion_damping_curve.Sgap_over_S_table;
        j["cushion_damping_curve"]["D_table"] = v.cushion_damping_curve.D_table;
    } else {
        j["D_damping"] = v.D_damping;
    }

    // Computed equilibrium — informational
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

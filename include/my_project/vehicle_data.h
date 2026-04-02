#pragma once
#include <string>

namespace acv {

struct VehicleData {
    // Constants
    double g = 9.81;
    double rho_air = 1.225;
    double p_atm = 101325.0;
    double n = 1.4;
    double chi = 0.6;

    // Geometry (half-vehicle)
    double m = 2050.0;      // kg
    double L = 10.5;        // m
    double l = 0.7;         // m
    double S = 20.0;        // m²
    double W0 = 14.0;       // m³
    double I_phi = 9.35;    // kg·m²

    // Equilibrium
    double p0 = 1005.5;     // Pa
    double Q0 = 5.71;       // m³/s
    double phi0 = 0.19;     // rad

    // Derivatives
    double dQin_dp = -0.002;
    double dM_dphi_per_L = -2178.0;  // N·m/(rad·m)
    double dM_dphidot = -3.0;        // N·m·s/rad
    double dY_dHdot_factor = -17836.0;

    // Helpers
    double get_dM_dphi() const { return dM_dphi_per_L * L; }
    double get_dY_dHdot() const { return dY_dHdot_factor * Q0 * rho_air; }
};

VehicleData loadVehicleFromJson(const std::string& filename);

} // namespace acv

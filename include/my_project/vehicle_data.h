#pragma once
#include <string>

namespace acv {

struct VehicleData {
    std::string name;

    // ===== Geometry =====
    double m;      // kg
    double L;      // m
    double l;      // m
    double S;      // m²
    double W0;     // m³
    double I_phi;  // kg·m²

    // ===== Equilibrium =====
    double p0;    // Pa
    double Q0;    // m³/s
    double phi0;  // rad

    // ===== Derivatives input =====
    double dQin_dp;
    double dM_dphi_per_L;
    double dM_dphidot;
    double dY_dHdot_factor;

    // ===== Derivatives computable =====
    double dQout_dp;
    double dQout_dphi;
    double dM_dphi;
    double dY_dHdot;

    // ===== API =====
    void Init();

   private:
    void Validate() const;
    void ComputeDerivatives();
};

VehicleData loadVehicleFromJson(const std::string& filename);

}  // namespace acv

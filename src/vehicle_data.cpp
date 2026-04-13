#include <cassert>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>

#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {

void VehicleData::Init() {
    Validate();
    ComputeDerivatives();
}

void VehicleData::Validate() const {
    // ===== Geometry =====
    assert(m > 0);
    assert(L > 0);
    assert(l > 0);
    assert(S > 0);
    assert(W0 > 0);
    assert(I_phi > 0);

    // ===== Equilibrium =====
    assert(p0 > 0);
    assert(Q0 > 0);
    assert(std::abs(phi0) < 1.5);

    // ===== Derivatives =====
    assert(std::isfinite(dQin_dp));
    assert(dQin_dp < 0);

    assert(std::isfinite(dM_dphi_per_L));
    assert(dM_dphi_per_L < 0);

    assert(std::isfinite(dM_dphidot));
    assert(dM_dphidot < 0);

    assert(std::isfinite(dY_dHdot_factor));
    assert(dY_dHdot_factor < 0);
}

void VehicleData::ComputeDerivatives() {
    using namespace constants;

    // ===== Outflow =====
    dQout_dp = Q0 / (2.0 * p0);

    dQout_dphi = kChi * std::sqrt(2.0 * p0 / kRhoAir) * L * l * std::sin(phi0);

    // ===== Seal =====
    dM_dphi = dM_dphi_per_L * L;

    // ===== Cushion =====
    dY_dHdot = dY_dHdot_factor;

    // ===== sanity =====
    assert(std::isfinite(dQout_dp));
    assert(std::isfinite(dQout_dphi));

    assert(std::isfinite(dM_dphi));
    assert(dM_dphi < 0);

    assert(std::isfinite(dY_dHdot));
    assert(dY_dHdot < 0);

    // ===== устойчивость (важно) =====
    double stability_term = dQin_dp - Q0 / (2.0 * p0);
    assert(stability_term < 0);
}

VehicleData loadVehicleFromJson(const std::string& filename) {
    std::ifstream file(filename);
    nlohmann::json j;
    file >> j;

    VehicleData v;

    v.name = j["name"];

    // ===== Geometry =====
    v.m = j["geometry"]["m"];
    v.L = j["geometry"]["L"];
    v.l = j["geometry"]["l"];
    v.S = j["geometry"]["S"];
    v.W0 = j["geometry"]["W0"];
    v.I_phi = j["geometry"]["I_phi"];

    // ===== Equilibrium =====
    v.p0 = j["equilibrium"]["p0"];
    v.Q0 = j["equilibrium"]["Q0"];
    v.phi0 = j["equilibrium"]["phi0"];

    // ===== Derivatives =====
    v.dQin_dp = j["derivatives"]["dQin_dp"];
    v.dM_dphi_per_L = j["derivatives"]["dM_dphi_per_L"];
    v.dM_dphidot = j["derivatives"]["dM_dphidot"];
    v.dY_dHdot_factor = j["derivatives"]["dY_dHdot_factor"];

    v.Init();

    return v;
}

}  // namespace acv

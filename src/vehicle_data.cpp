#include "my_project/vehicle_data.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

using json = nlohmann::json;

namespace acv {

VehicleData loadVehicleFromJson(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    json j;
    file >> j;

    VehicleData v;

    // Constants
    v.g = j["constants"]["g"];
    v.rho_air = j["constants"]["rho_air"];
    v.p_atm = j["constants"]["p_atm"];
    v.n = j["constants"]["n_polytropic"];
    v.chi = j["constants"]["chi_outflow"];

    // Geometry
    v.m = j["geometry"]["m_half_kg"];
    v.L = j["geometry"]["L_m"];
    v.l = j["geometry"]["l_m"];
    v.S = j["geometry"]["S_m2"];
    v.W0 = j["geometry"]["W0_m3"];
    v.I_phi = j["geometry"]["I_phi_kgm2"];

    // Equilibrium
    v.p0 = j["equilibrium"]["p0_Pa"];
    v.Q0 = j["equilibrium"]["Q0_m3s"];
    v.phi0 = j["equilibrium"]["phi0_rad"];

    // Derivatives
    v.dQin_dp = j["derivatives"]["dQin_dp"];
    v.dM_dphi_per_L = j["derivatives"]["dM_dphi_per_L"];
    v.dM_dphidot = j["derivatives"]["dM_dphidot"];
    v.dY_dHdot_factor = j["derivatives"]["dY_dHdot_factor"];

    return v;
}

} // namespace acv

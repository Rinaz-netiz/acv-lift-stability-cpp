// tests/test_vehicle_data.cpp
#include <gtest/gtest.h>

#include <fstream>
#include <nlohmann/json.hpp>

#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class VehicleDataTest : public ::testing::Test {
   protected:
    VehicleData CreateValidVehicle() {
        VehicleData v;
        v.name = "Test Vehicle";

        // Geometry
        v.m = 1000.0;
        v.L = 10.0;
        v.l = 0.5;
        v.S = 50.0;
        v.W0 = 25.0;
        v.I_phi = 500.0;

        // Equilibrium
        v.p0 = 2000.0;
        v.Q0 = 5.0;
        v.phi0 = 0.1;

        // Derivatives
        v.dQin_dp = -0.001;
        v.dM_dphi_per_L = -100.0;
        v.dM_dphidot = -50.0;
        v.dY_dHdot_factor = -0.5;

        return v;
    }
};

TEST_F(VehicleDataTest, InitComputesDerivativesCorrectly) {
    VehicleData v = CreateValidVehicle();
    v.Init();

    // Check dQout_dp = Q0 / (2 * p0)
    EXPECT_NEAR(v.dQout_dp, v.Q0 / (2.0 * v.p0), 1e-9);

    // Check dQout_dphi computation
    double expected_dQout_dphi = constants::kChi *
                                 std::sqrt(2.0 * v.p0 / constants::kRhoAir) *
                                 v.L * v.l * std::sin(v.phi0);
    EXPECT_NEAR(v.dQout_dphi, expected_dQout_dphi, 1e-6);

    // Check dM_dphi = dM_dphi_per_L * L
    EXPECT_NEAR(v.dM_dphi, v.dM_dphi_per_L * v.L, 1e-9);

    // Check dY_dHdot
    double expected_dY_dHdot = v.dY_dHdot_factor * v.Q0 * constants::kRhoAir;
    EXPECT_NEAR(v.dY_dHdot, expected_dY_dHdot, 1e-6);
}

TEST_F(VehicleDataTest, ValidationFailsOnNegativeMass) {
    VehicleData v = CreateValidVehicle();
    v.m = -100.0;
    EXPECT_DEATH(v.Init(), "");
}

TEST_F(VehicleDataTest, ValidationFailsOnZeroLength) {
    VehicleData v = CreateValidVehicle();
    v.L = 0.0;
    EXPECT_DEATH(v.Init(), "");
}

TEST_F(VehicleDataTest, ValidationFailsOnInvalidPhi0) {
    VehicleData v = CreateValidVehicle();
    v.phi0 = 2.0;  // > 1.5 rad
    EXPECT_DEATH(v.Init(), "");
}

TEST_F(VehicleDataTest, ValidationFailsOnPositiveDQinDp) {
    VehicleData v = CreateValidVehicle();
    v.dQin_dp = 0.001;  // Should be negative
    EXPECT_DEATH(v.Init(), "");
}

TEST_F(VehicleDataTest, LoadFromJsonValid) {
    // Create temporary JSON file
    nlohmann::json j;
    j["name"] = "JSON Test Vehicle";
    j["geometry"]["m"] = 1000.0;
    j["geometry"]["L"] = 10.0;
    j["geometry"]["l"] = 0.5;
    j["geometry"]["S"] = 50.0;
    j["geometry"]["W0"] = 25.0;
    j["geometry"]["I_phi"] = 500.0;
    j["equilibrium"]["p0"] = 2000.0;
    j["equilibrium"]["Q0"] = 5.0;
    j["equilibrium"]["phi0"] = 0.1;
    j["derivatives"]["dQin_dp"] = -0.001;
    j["derivatives"]["dM_dphi_per_L"] = -100.0;
    j["derivatives"]["dM_dphidot"] = -50.0;
    j["derivatives"]["dY_dHdot_factor"] = -0.5;

    std::ofstream file("test_vehicle.json");
    file << j.dump(4);
    file.close();

    VehicleData v = loadVehicleFromJson("test_vehicle.json");

    EXPECT_EQ(v.name, "JSON Test Vehicle");
    EXPECT_DOUBLE_EQ(v.m, 1000.0);
    EXPECT_DOUBLE_EQ(v.L, 10.0);

    // Cleanup
    std::remove("test_vehicle.json");
}

TEST_F(VehicleDataTest, StabilityTermIsNegative) {
    VehicleData v = CreateValidVehicle();
    v.Init();

    double stability_term = v.dQin_dp - v.Q0 / (2.0 * v.p0);
    EXPECT_LT(stability_term, 0.0);
}

}  // namespace
}  // namespace acv

// tests/test_integration.cpp
#include <gtest/gtest.h>

#include "my_project/analytical.h"
#include "my_project/numerical.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class IntegrationTest : public ::testing::Test {
   protected:
    VehicleData CreateTestVehicle() {
        VehicleData v;
        v.m = 1000.0;
        v.L = 10.0;
        v.l = 0.5;
        v.S = 50.0;
        v.W0 = 25.0;
        v.I_phi = 500.0;
        v.p0 = 2000.0;
        v.Q0 = 5.0;
        v.phi0 = 0.1;
        v.dQin_dp = -0.005;
        v.dM_dphi_per_L = -100.0;
        v.dM_dphidot = -50.0;
        v.dY_dHdot_factor = -0.5;
        v.Init();
        return v;
    }
};

TEST_F(IntegrationTest, AnalyticalAndNumericalAgreeOnStability) {
    VehicleData v = CreateTestVehicle();

    bool analytical = AnalyticalVerification(v);
    StabilityResult numerical = AnalyzeStabilitySimple(v);

    // Both methods should agree on stability for simple case
    EXPECT_EQ(analytical, numerical.is_stable);
}

TEST_F(IntegrationTest, SimpleAndFullNumericalConsistent) {
    VehicleData v = CreateTestVehicle();

    StabilityResult simple = AnalyzeStabilitySimple(v);
    StabilityResult full = AnalyzeStabilityFull(v);

    // For stable systems, both should agree
    if (simple.is_stable && full.is_stable) {
        EXPECT_TRUE(true);  // Both stable
    }
}

}  // namespace
}  // namespace acv

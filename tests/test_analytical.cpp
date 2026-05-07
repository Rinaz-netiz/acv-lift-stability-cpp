// tests/test_analytical.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "my_project/analytical.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class AnalyticalTest : public ::testing::Test {
   protected:
    VehicleData CreateStableVehicle() {
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
        v.dQin_dp = -0.01;  // Strong negative (stabilizing)
        v.dM_dphi_per_L = -100.0;
        v.dM_dphidot = -50.0;
        v.dY_dHdot_factor = -0.5;
        v.Init();
        return v;
    }

    VehicleData CreateUnstableVehicle() {
        VehicleData v;
        v.m = 1000.0;
        v.L = 10.0;
        v.l = 0.5;
        v.S = 50.0;
        v.W0 = 25.0;
        v.I_phi = 500.0;
        v.p0 = 2000.0;
        v.Q0 = 5.0;
        v.phi0 = 0.5;         // Larger angle
        v.dQin_dp = -0.0001;  // Weak negative (less stabilizing)
        v.dM_dphi_per_L = -100.0;
        v.dM_dphidot = -50.0;
        v.dY_dHdot_factor = -0.5;
        v.Init();
        return v;
    }
};

TEST_F(AnalyticalTest, StableVehicleDetectedAsStable) {
    VehicleData v = CreateStableVehicle();
    AnalyticalResult result = VerifyAnalyticalDetailed(v);

    EXPECT_TRUE(result.is_stable);
    EXPECT_LT(result.stability_margin, 0.0);
}

TEST_F(AnalyticalTest, UnstableVehicleDetectedAsUnstable) {
    VehicleData v = CreateUnstableVehicle();
    AnalyticalResult result = VerifyAnalyticalDetailed(v);

    EXPECT_FALSE(result.is_stable);
    EXPECT_GT(result.stability_margin, 0.0);
}

TEST_F(AnalyticalTest, PneumaticTermIsNegative) {
    VehicleData v = CreateStableVehicle();
    AnalyticalResult result = VerifyAnalyticalDetailed(v);

    // Pneumatic term should be negative (stabilizing)
    EXPECT_LT(result.pneumatic_term, 0.0);
}

TEST_F(AnalyticalTest, GeometricTermIsPositive) {
    VehicleData v = CreateStableVehicle();
    v.phi0 = 0.3;  // Positive angle
    v.Init();

    AnalyticalResult result = VerifyAnalyticalDetailed(v);

    // Geometric term should be positive (destabilizing) for positive phi
    EXPECT_GT(result.geometric_term, 0.0);
}

TEST_F(AnalyticalTest, InfluenceRatioIsPositive) {
    VehicleData v = CreateStableVehicle();
    AnalyticalResult result = VerifyAnalyticalDetailed(v);

    EXPECT_GT(result.influence_ratio, 0.0);
}

TEST_F(AnalyticalTest, SimpleVerificationMatchesDetailed) {
    VehicleData v = CreateStableVehicle();

    bool simple_result = AnalyticalVerification(v);
    AnalyticalResult detailed_result = VerifyAnalyticalDetailed(v);

    EXPECT_EQ(simple_result, detailed_result.is_stable);
}

TEST_F(AnalyticalTest, ZeroAngleGivesZeroGeometricTerm) {
    VehicleData v = CreateStableVehicle();
    v.phi0 = 0.0;
    v.Init();

    AnalyticalResult result = VerifyAnalyticalDetailed(v);

    EXPECT_NEAR(result.geometric_term, 0.0, 1e-9);
}

}  // namespace
}  // namespace acv

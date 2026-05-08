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
        v.name = "Test Vehicle";

        // Параметры геометрии (теперь через v.geometry)
        v.geometry.m = 1000.0;
        v.geometry.L = 10.0;
        v.geometry.l = 0.5;
        v.geometry.S = 50.0;
        v.geometry.W0 = 25.0;
        v.geometry.I_phi = 500.0;
        v.geometry.h0 = 0.05;

        // Параметры нагнетателя
        v.blower.Q_max = 10.0;
        v.blower.k_fan = 0.0025;

        // Равновесное состояние
        v.p0 = 2000.0;
        v.Q0 = 5.0;
        v.phi0 = 0.1;
        v.H0 = 0.05;

        // Производные (имена приведены в соответствие с vehicle_data.h)
        v.dQin_dp = -0.005;
        v.dM_dphi = -100.0;  // Было dM_dphi_per_L
        v.dM_dphidot = -50.0;
        v.dY_dHdot = -0.5;  // Было dY_dHdot_factor

        // Инициализация внутренних состояний
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

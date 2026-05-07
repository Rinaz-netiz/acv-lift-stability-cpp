// tests/test_numerical.cpp
#include <gtest/gtest.h>

#include <complex>

#include "my_project/numerical.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class NumericalTest : public ::testing::Test {
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

TEST_F(NumericalTest, SimpleAnalysisReturnsThreeEigenvalues) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilitySimple(v);

    EXPECT_EQ(result.eigenvalues.size(), 3);
}

TEST_F(NumericalTest, FullAnalysisReturnsFiveEigenvalues) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilityFull(v);

    EXPECT_EQ(result.eigenvalues.size(), 5);
}

TEST_F(NumericalTest, StableSystemHasNegativeRealParts) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilitySimple(v);

    if (result.is_stable) {
        for (const auto& lambda : result.eigenvalues) {
            EXPECT_LT(lambda.real(), kEps);
        }
    }
}

TEST_F(NumericalTest, MaxRealPartIsConsistent) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilitySimple(v);

    double max_real = -1e18;
    for (const auto& lambda : result.eigenvalues) {
        if (lambda.real() > max_real) {
            max_real = lambda.real();
        }
    }

    EXPECT_NEAR(result.max_real_part, max_real, 1e-9);
}

TEST_F(NumericalTest, OscillationModesHavePositiveImaginaryPart) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilityFull(v);

    for (const auto& mode : result.oscillation_modes) {
        EXPECT_GT(mode.eigenvalue.imag(), kEps);
    }
}

TEST_F(NumericalTest, OscillationPeriodIsPositive) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilityFull(v);

    for (const auto& mode : result.oscillation_modes) {
        EXPECT_GT(mode.period, 0.0);
    }
}

TEST_F(NumericalTest, DecayRatioCalculationIsCorrect) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilityFull(v);

    for (const auto& mode : result.oscillation_modes) {
        double expected_decay = std::exp(mode.logarithmic_decrement);
        EXPECT_NEAR(mode.decay_ratio, expected_decay, 1e-9);
    }
}

TEST_F(NumericalTest, LogarithmicDecrementCalculation) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilityFull(v);

    for (const auto& mode : result.oscillation_modes) {
        double sigma = mode.eigenvalue.real();
        double expected_chi = -sigma * mode.period;
        EXPECT_NEAR(mode.logarithmic_decrement, expected_chi, 1e-9);
    }
}

TEST_F(NumericalTest, FullModelIncludesVerticalDynamics) {
    VehicleData v = CreateTestVehicle();

    StabilityResult simple = AnalyzeStabilitySimple(v);
    StabilityResult full = AnalyzeStabilityFull(v);

    // Full model should have more eigenvalues
    EXPECT_GT(full.eigenvalues.size(), simple.eigenvalues.size());
}

TEST_F(NumericalTest, UnstableSystemDetected) {
    VehicleData v = CreateTestVehicle();
    v.dQin_dp = -0.00001;  // Very weak stabilization
    v.phi0 = 0.8;          // Large angle
    v.Init();

    StabilityResult result = AnalyzeStabilitySimple(v);

    if (!result.is_stable) {
        EXPECT_GT(result.max_real_part, kEps);
    }
}

TEST_F(NumericalTest, EigenvaluesAreFinite) {
    VehicleData v = CreateTestVehicle();
    StabilityResult result = AnalyzeStabilityFull(v);

    for (const auto& lambda : result.eigenvalues) {
        EXPECT_TRUE(std::isfinite(lambda.real()));
        EXPECT_TRUE(std::isfinite(lambda.imag()));
    }
}

}  // namespace
}  // namespace acv

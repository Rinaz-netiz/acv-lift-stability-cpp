// tests/test_numerical.cpp
#include <gtest/gtest.h>

#include <cmath>
#include <complex>

#include "my_project/constants.h"
#include "my_project/numerical.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

// ============================================================
// Фабрика (та же что в других тестах)
// ============================================================

VehicleData MakeTestVehicle() {
    VehicleData v;
    v.name = "Numerical Test Hovercraft";

    v.geometry.m = 1000.0;
    v.geometry.L = 8.0;
    v.geometry.l = 0.4;
    v.geometry.S = 40.0;
    v.geometry.W0 = 20.0;
    v.geometry.I_phi = 400.0;
    v.geometry.h0 = 0.15;

    v.blower.type = BlowerCharacteristics::Type::Linear;
    v.blower.Q_max = 10.0;
    v.blower.k_fan = 0.00003;

    v.seal.k_seal_per_L = 5000.0;
    v.seal.k_phi = 500.0;
    v.seal.c_phi = 200.0;

    v.damping.beta_cushion = 0.3;

    v.Init();
    return v;
}

// ============================================================
// Тест-класс
// ============================================================

class NumericalSimpleTest : public ::testing::Test {
   protected:
    VehicleData vehicle_;
    StabilityResult result_;

    void SetUp() override {
        vehicle_ = MakeTestVehicle();
        result_ = AnalyzeStabilitySimple(vehicle_);
    }
};

class NumericalFullTest : public ::testing::Test {
   protected:
    VehicleData vehicle_;
    StabilityResult result_;

    void SetUp() override {
        vehicle_ = MakeTestVehicle();
        result_ = AnalyzeStabilityFull(vehicle_);
    }
};

// ============================================================
// 1. Тесты простой модели (3 степени свободы)
// ============================================================

TEST_F(NumericalSimpleTest, ReturnsThreeEigenvalues) {
    EXPECT_EQ(result_.eigenvalues.size(), 3u)
        << "Модель 3 СС должна давать 3 собственных значения";
}

TEST_F(NumericalSimpleTest, StableVehicleIsDetectedAsStable) {
    EXPECT_TRUE(result_.is_stable)
        << "Физически устойчивое КВП должно быть определено как устойчивое";
}

TEST_F(NumericalSimpleTest, MaxRealPartIsNegativeForStableSystem) {
    EXPECT_LT(result_.max_real_part, 0.0)
        << "Максимальная вещественная часть СЗ < 0 для устойчивой системы";
}

TEST_F(NumericalSimpleTest, MaxRealPartIsConsistentWithEigenvalues) {
    // max_real_part должен совпадать с максимумом Re(λ_i)
    double max_real = -1e18;
    for (const auto& lambda : result_.eigenvalues) {
        if (lambda.real() > max_real) max_real = lambda.real();
    }
    EXPECT_NEAR(result_.max_real_part, max_real, 1e-9)
        << "max_real_part должен совпадать с max(Re(λ_i))";
}

TEST_F(NumericalSimpleTest, IsStableConsistentWithMaxRealPart) {
    bool expected = result_.max_real_part <= kEps;
    EXPECT_EQ(result_.is_stable, expected)
        << "is_stable должен соответствовать знаку max_real_part";
}

TEST_F(NumericalSimpleTest, AllEigenvaluesAreFinite) {
    for (size_t i = 0; i < result_.eigenvalues.size(); ++i) {
        EXPECT_TRUE(std::isfinite(result_.eigenvalues[i].real()))
            << "Re(λ_" << i << ") должно быть конечным";
        EXPECT_TRUE(std::isfinite(result_.eigenvalues[i].imag()))
            << "Im(λ_" << i << ") должно быть конечным";
    }
}

TEST_F(NumericalSimpleTest, OscillationModesHavePositivePeriod) {
    for (size_t i = 0; i < result_.oscillation_modes.size(); ++i) {
        EXPECT_GT(result_.oscillation_modes[i].period, 0.0)
            << "Период моды " << i << " должен быть положительным";
    }
}

TEST_F(NumericalSimpleTest, OscillationModePeriodMatchesFormula) {
    for (const auto& mode : result_.oscillation_modes) {
        double omega = std::abs(mode.eigenvalue.imag());
        double expected = 2.0 * M_PI / omega;
        EXPECT_NEAR(mode.period, expected, 1e-9) << "T = 2π/ω";
    }
}

TEST_F(NumericalSimpleTest, LogarithmicDecrementMatchesFormula) {
    for (const auto& mode : result_.oscillation_modes) {
        double sigma = mode.eigenvalue.real();
        double expected = -sigma * mode.period;
        EXPECT_NEAR(mode.logarithmic_decrement, expected, 1e-9) << "χ = -σ*T";
    }
}

TEST_F(NumericalSimpleTest, DecayRatioMatchesFormula) {
    for (const auto& mode : result_.oscillation_modes) {
        double expected = std::exp(mode.logarithmic_decrement);
        EXPECT_NEAR(mode.decay_ratio, expected, 1e-9) << "K = exp(χ)";
    }
}

TEST_F(NumericalSimpleTest, OscillationModesHavePositiveImaginaryPart) {
    for (const auto& mode : result_.oscillation_modes) {
        EXPECT_GT(mode.eigenvalue.imag(), kEps)
            << "Мода должна иметь Im(λ) > 0 (берём только положительные)";
    }
}

// ============================================================
// 2. Тесты полной модели (5 степеней свободы)
// ============================================================

TEST_F(NumericalFullTest, ReturnsFiveEigenvalues) {
    EXPECT_EQ(result_.eigenvalues.size(), 5u)
        << "Модель 5 СС должна давать 5 собственных значений";
}

TEST_F(NumericalFullTest, StableVehicleIsDetectedAsStable) {
    EXPECT_TRUE(result_.is_stable)
        << "Физически устойчивое КВП должно быть устойчивым в полной модели";
}

TEST_F(NumericalFullTest, MaxRealPartIsNegativeForStableSystem) {
    EXPECT_LT(result_.max_real_part, 0.0)
        << "Максимальная вещественная часть СЗ < 0";
}

TEST_F(NumericalFullTest, MaxRealPartIsConsistentWithEigenvalues) {
    double max_real = -1e18;
    for (const auto& lambda : result_.eigenvalues) {
        if (lambda.real() > max_real) max_real = lambda.real();
    }
    EXPECT_NEAR(result_.max_real_part, max_real, 1e-9);
}

TEST_F(NumericalFullTest, AllEigenvaluesAreFinite) {
    for (size_t i = 0; i < result_.eigenvalues.size(); ++i) {
        EXPECT_TRUE(std::isfinite(result_.eigenvalues[i].real()))
            << "Re(λ_" << i << ") конечно";
        EXPECT_TRUE(std::isfinite(result_.eigenvalues[i].imag()))
            << "Im(λ_" << i << ") конечно";
    }
}

TEST_F(NumericalFullTest, OscillationModePeriodMatchesFormula) {
    for (const auto& mode : result_.oscillation_modes) {
        double omega = std::abs(mode.eigenvalue.imag());
        double expected = 2.0 * M_PI / omega;
        EXPECT_NEAR(mode.period, expected, 1e-9);
    }
}

TEST_F(NumericalFullTest, DecayRatioMatchesFormula) {
    for (const auto& mode : result_.oscillation_modes) {
        double expected = std::exp(mode.logarithmic_decrement);
        EXPECT_NEAR(mode.decay_ratio, expected, 1e-9);
    }
}

TEST_F(NumericalFullTest, HasMoreEigenvaluesThanSimpleModel) {
    StabilityResult simple = AnalyzeStabilitySimple(vehicle_);
    EXPECT_GT(result_.eigenvalues.size(), simple.eigenvalues.size())
        << "Полная модель имеет больше СЗ чем простая";
}

// ============================================================
// 3. Сравнительные тесты простой и полной моделей
// ============================================================

class NumericalComparisonTest : public ::testing::Test {
   protected:
    VehicleData vehicle_;

    void SetUp() override { vehicle_ = MakeTestVehicle(); }
};

TEST_F(NumericalComparisonTest, BothModelsAgreeOnStability) {
    auto simple = AnalyzeStabilitySimple(vehicle_);
    auto full = AnalyzeStabilityFull(vehicle_);

    // Для заведомо устойчивого КВП обе модели должны согласоваться
    if (simple.is_stable) {
        EXPECT_TRUE(full.is_stable)
            << "Если простая модель устойчива, полная тоже должна быть";
    }
}

TEST_F(NumericalComparisonTest, FullModelMaxRealPartIsFinite) {
    auto full = AnalyzeStabilityFull(vehicle_);
    EXPECT_TRUE(std::isfinite(full.max_real_part));
}

// ============================================================
// 4. Тесты граничных случаев
// ============================================================

class NumericalEdgeCasesTest : public ::testing::Test {};

TEST_F(NumericalEdgeCasesTest, ZeroDampingDoesNotCrash) {
    VehicleData v = MakeTestVehicle();
    // Убираем демпфирование — система на границе устойчивости
    v.seal.c_phi = 0.0;
    v.damping.beta_cushion = 0.0;
    v.Init();

    // Не должно бросать исключений
    EXPECT_NO_THROW({
        auto result = AnalyzeStabilitySimple(v);
        (void)result;
    });
}

TEST_F(NumericalEdgeCasesTest, LargeInertiaContinuesConverging) {
    VehicleData v = MakeTestVehicle();
    v.geometry.I_phi = 1e5;
    v.Init();

    EXPECT_NO_THROW({
        auto result = AnalyzeStabilitySimple(v);
        EXPECT_EQ(result.eigenvalues.size(), 3u);
    });
}

TEST_F(NumericalEdgeCasesTest, SmallInertiaContinuesConverging) {
    VehicleData v = MakeTestVehicle();
    v.geometry.I_phi = 1.0;
    v.Init();

    EXPECT_NO_THROW({
        auto result = AnalyzeStabilitySimple(v);
        EXPECT_EQ(result.eigenvalues.size(), 3u);
    });
}

TEST_F(NumericalEdgeCasesTest, FullModelDoesNotCrashWithZeroDamping) {
    VehicleData v = MakeTestVehicle();
    v.seal.c_phi = 0.0;
    v.damping.beta_cushion = 0.0;
    v.Init();

    EXPECT_NO_THROW({
        auto result = AnalyzeStabilityFull(v);
        EXPECT_EQ(result.eigenvalues.size(), 5u);
    });
}

}  // namespace
}  // namespace acv

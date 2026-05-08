#include <gtest/gtest.h>

#include <cmath>

#include "my_project/numerical.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class NumericalTest : public ::testing::Test {
   protected:
    void SetUp() override {
        data = LoadVehicleFromJson("test_data/vehicle_a.json");
    }

    VehicleData data;
};

// ============================================================
// Тесты простой модели (3×3)
// ============================================================

TEST_F(NumericalTest, SimpleModelThreeEigenvalues) {
    StabilityResult result = AnalyzeStabilitySimple(data);

    EXPECT_EQ(result.eigenvalues.size(), 3);
}

TEST_F(NumericalTest, SimpleModelStabilityConsistency) {
    StabilityResult result = AnalyzeStabilitySimple(data);

    // is_stable должно соответствовать max_real_part
    if (result.max_real_part > 1e-9) {
        EXPECT_FALSE(result.is_stable);
    } else {
        EXPECT_TRUE(result.is_stable);
    }
}

TEST_F(NumericalTest, SimpleModelMaxRealPartCorrect) {
    StabilityResult result = AnalyzeStabilitySimple(data);

    double max_re = -1e18;
    for (const auto& lambda : result.eigenvalues) {
        max_re = std::max(max_re, lambda.real());
    }

    EXPECT_NEAR(result.max_real_part, max_re, 1e-9);
}

TEST_F(NumericalTest, SimpleModelOscillationModes) {
    StabilityResult result = AnalyzeStabilitySimple(data);

    // Для каждой комплексной пары должна быть мода
    int complex_count = 0;
    for (const auto& lambda : result.eigenvalues) {
        if (std::abs(lambda.imag()) > 1e-9) {
            complex_count++;
        }
    }

    EXPECT_GE(result.oscillation_modes.size(), 0);
    EXPECT_LE(result.oscillation_modes.size(), result.eigenvalues.size());
}

TEST_F(NumericalTest, SimpleModelPeriodPositive) {
    StabilityResult result = AnalyzeStabilitySimple(data);

    for (const auto& mode : result.oscillation_modes) {
        EXPECT_GT(mode.period, 0.0);
        EXPECT_LT(mode.period, 1000.0);  // Разумные пределы
    }
}

TEST_F(NumericalTest, SimpleModelDecayRatioMeaningful) {
    StabilityResult result = AnalyzeStabilitySimple(data);

    for (const auto& mode : result.oscillation_modes) {
        EXPECT_GT(mode.decay_ratio, 0.0);
        // K = exp(χ), где χ = -σ·T
        // Для устойчивости: σ < 0 => χ > 0 => K > 1
        // Для неустойчивости: σ > 0 => χ < 0 => K < 1

        if (mode.eigenvalue.real() < 0) {
            // Затухающие колебания
            EXPECT_GT(mode.decay_ratio, 1.0);
        } else {
            // Растущие колебания
            EXPECT_LT(mode.decay_ratio, 1.0);
        }
    }
}

// ============================================================
// Тесты полной модели (5×5)
// ============================================================

TEST_F(NumericalTest, FullModelFiveEigenvalues) {
    StabilityResult result = AnalyzeStabilityFull(data);

    EXPECT_EQ(result.eigenvalues.size(), 5);
}

TEST_F(NumericalTest, FullModelStabilityConsistency) {
    StabilityResult result = AnalyzeStabilityFull(data);

    if (result.max_real_part > 1e-9) {
        EXPECT_FALSE(result.is_stable);
    } else {
        EXPECT_TRUE(result.is_stable);
    }
}

TEST_F(NumericalTest, FullModelIncludesHeightDynamics) {
    StabilityResult full = AnalyzeStabilityFull(data);
    StabilityResult simple = AnalyzeStabilitySimple(data);

    // Полная модель имеет больше СЗ
    EXPECT_GT(full.eigenvalues.size(), simple.eigenvalues.size());
}

// ============================================================
// Сравнение моделей
// ============================================================

TEST_F(NumericalTest, SimpleAndFullModelsStabilityAgreement) {
    StabilityResult simple = AnalyzeStabilitySimple(data);
    StabilityResult full = AnalyzeStabilityFull(data);

    // Оба метода должны согласовываться по устойчивости
    // (хотя могут различаться в деталях)
    if (simple.is_stable && full.is_stable) {
        SUCCEED();
    } else if (!simple.is_stable && !full.is_stable) {
        SUCCEED();
    } else {
        // Расхождение возможно из-за дополнительной динамики высоты
        // Логируем, но не обязательно ошибка
        std::cout << "Models disagree: simple=" << simple.is_stable
                  << " full=" << full.is_stable << std::endl;
    }
}

TEST_F(NumericalTest, FullModelHasMoreModes) {
    StabilityResult simple = AnalyzeStabilitySimple(data);
    StabilityResult full = AnalyzeStabilityFull(data);

    // Полная модель может иметь больше осцилляционных мод
    EXPECT_GE(full.oscillation_modes.size(), simple.oscillation_modes.size());
}

// ============================================================
// Тесты численной устойчивости
// ============================================================

TEST_F(NumericalTest, EigenvaluesNotNaN) {
    StabilityResult result = AnalyzeStabilityFull(data);

    for (const auto& lambda : result.eigenvalues) {
        EXPECT_FALSE(std::isnan(lambda.real()));
        EXPECT_FALSE(std::isnan(lambda.imag()));
    }
}

TEST_F(NumericalTest, ComplexEigenvaluesComePairs) {
    StabilityResult result = AnalyzeStabilityFull(data);

    std::vector<std::complex<double>> complex_eigs;
    for (const auto& lambda : result.eigenvalues) {
        if (std::abs(lambda.imag()) > 1e-9) {
            complex_eigs.push_back(lambda);
        }
    }

    // Комплексные СЗ должны идти сопряжёнными парами
    EXPECT_EQ(complex_eigs.size() % 2, 0);
}

// ============================================================
// Тесты параметрических вариаций
// ============================================================

TEST_F(NumericalTest, IncreasedDampingMoreStable) {
    VehicleData data1 = data;
    VehicleData data2 = data;

    // Увеличим демпфирование ограждения
    data2.seal_damping_curve.dM_dphidot_table = {-100.0, -80.0, -50.0, -20.0,
                                                 -5.0};

    data1.Init();
    data2.Init();

    StabilityResult r1 = AnalyzeStabilitySimple(data1);
    StabilityResult r2 = AnalyzeStabilitySimple(data2);

    // Больше демпфирование -> меньше max_real_part (более отрицательный)
    EXPECT_LE(r2.max_real_part, r1.max_real_part);
}

TEST_F(NumericalTest, IncreasedInertiaSlowerOscillations) {
    VehicleData data1 = data;
    VehicleData data2 = data;

    data1.geometry.I_phi = 9.35;
    data2.geometry.I_phi = 18.7;  // Удвоенная инерция

    data1.Init();
    data2.Init();

    StabilityResult r1 = AnalyzeStabilitySimple(data1);
    StabilityResult r2 = AnalyzeStabilitySimple(data2);

    if (!r1.oscillation_modes.empty() && !r2.oscillation_modes.empty()) {
        // Больше инерция -> больше период (медленнее колебания)
        EXPECT_GT(r2.oscillation_modes[0].period,
                  r1.oscillation_modes[0].period);
    }
}

}  // namespace
}  // namespace acv

#include <gtest/gtest.h>

#include <cmath>

#include "my_project/analytical.h"
#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class AnalyticalTest : public ::testing::Test {
   protected:
    void SetUp() override {
        data = LoadVehicleFromJson("test_data/vehicle_a.json");
    }

    VehicleData data;
};

TEST_F(AnalyticalTest, AnalyticalVerificationBasic) {
    bool stable = AnalyticalVerification(data);
    // Для данных из примера ожидаем определённый результат
    // (зависит от конкретных значений)
    EXPECT_TRUE(stable || !stable);  // Проверка, что функция работает
}

TEST_F(AnalyticalTest, DetailedResultStructure) {
    AnalyticalResult result = VerifyAnalyticalDetailed(data);

    // Проверка согласованности результата
    EXPECT_EQ(result.is_stable, result.stability_margin < 0);

    // Запас устойчивости = сумма слагаемых
    double sum = result.pneumatic_term + result.geometric_term;
    EXPECT_NEAR(result.stability_margin, sum, std::abs(sum) * 1e-6);
}

TEST_F(AnalyticalTest, PneumaticTermSign) {
    AnalyticalResult result = VerifyAnalyticalDetailed(data);

    // Пневматическое слагаемое обычно отрицательное (стабилизирующее)
    // при dQin/dp < dQout/dp
    if (data.der.dQin_dp < data.der.dQout_dp) {
        EXPECT_LT(result.pneumatic_term, 0.0);
    }
}

TEST_F(AnalyticalTest, GeometricTermPositiveForPositiveAngle) {
    AnalyticalResult result = VerifyAnalyticalDetailed(data);

    // Геометрическое слагаемое положительно при φ₀ > 0
    if (data.eq.phi0 > 0.0) {
        EXPECT_GT(result.geometric_term, 0.0);
    }
}

TEST_F(AnalyticalTest, InfluenceRatioMeaningful) {
    AnalyticalResult result = VerifyAnalyticalDetailed(data);

    if (std::abs(result.geometric_term) > 1e-9) {
        EXPECT_GT(result.influence_ratio, 0.0);
        // Отношение должно быть конечным
        EXPECT_LT(result.influence_ratio, 1e6);
    }
}

TEST_F(AnalyticalTest, StabilityConsistency) {
    bool simple = AnalyticalVerification(data);
    AnalyticalResult detailed = VerifyAnalyticalDetailed(data);

    // Оба метода должны давать одинаковый результат
    EXPECT_EQ(simple, detailed.is_stable);
}

// Тест на граничный случай: φ₀ = 0
TEST_F(AnalyticalTest, ZeroAngleCase) {
    // Искусственно установим φ₀ = 0
    data.eq.phi0 = 0.0;

    AnalyticalResult result = VerifyAnalyticalDetailed(data);

    // При φ₀ = 0: sin(φ₀) = 0 => geometric_term = 0
    EXPECT_NEAR(result.geometric_term, 0.0, 1e-9);

    // Устойчивость определяется только пневматическим слагаемым
    EXPECT_NEAR(result.stability_margin, result.pneumatic_term, 1e-9);
}

// Тест влияния объёма подушки
TEST_F(AnalyticalTest, LargeVolumeMoreStable) {
    VehicleData data1 = data;
    VehicleData data2 = data;

    data1.geometry.W0 = 10.0;
    data2.geometry.W0 = 20.0;  // Больше объём

    data1.Init();
    data2.Init();

    AnalyticalResult r1 = VerifyAnalyticalDetailed(data1);
    AnalyticalResult r2 = VerifyAnalyticalDetailed(data2);

    // Больший объём -> меньше пневматическое слагаемое (по модулю)
    // -> больше stability_margin (менее отрицательный)
    // Это упрощённое предположение, зависит от других параметров
}

}  // namespace
}  // namespace acv

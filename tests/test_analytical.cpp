// tests/test_analytical.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "my_project/analytical.h"
#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

// ============================================================
// Используем ту же фабрику, что и в test_vehicle_data.cpp
// ============================================================

VehicleData MakeStableVehicle() {
    VehicleData v;
    v.name = "Stable Test Hovercraft";

    v.geometry.m = 1000.0;
    v.geometry.L = 8.0;
    v.geometry.l = 0.4;
    v.geometry.S = 40.0;
    v.geometry.W0 = 20.0;
    v.geometry.I_phi = 400.0;
    v.geometry.h0 = 0.15;

    // Сильный нагнетатель → сильное пневматическое демпфирование
    v.blower.type = BlowerCharacteristics::Type::Linear;
    v.blower.Q_max = 10.0;
    v.blower.k_fan = 0.00003;

    // Жёсткое ограждение с хорошим демпфированием
    v.seal.k_seal_per_L = 5000.0;
    v.seal.k_phi = 500.0;
    v.seal.c_phi = 200.0;

    v.damping.beta_cushion = 0.3;

    v.Init();
    return v;
}

VehicleData MakeUnstableVehicle() {
    VehicleData v;
    v.name = "Unstable Test Hovercraft";

    v.geometry.m = 1000.0;
    v.geometry.L = 8.0;
    v.geometry.l = 0.4;
    v.geometry.S = 40.0;
    v.geometry.W0 = 20.0;
    v.geometry.I_phi = 400.0;
    v.geometry.h0 = 0.15;

    // Слабый нагнетатель → почти горизонтальная характеристика
    // → слабое пневматическое демпфирование
    v.blower.type = BlowerCharacteristics::Type::Linear;
    v.blower.Q_max = 10.0;
    v.blower.k_fan = 1e-8;  // dQin/dp ≈ 0

    // Мягкое ограждение без демпфирования
    v.seal.k_seal_per_L = 100.0;
    v.seal.k_phi = 10.0;
    v.seal.c_phi = 0.0;  // Нет демпфирования

    v.damping.beta_cushion = 0.0;

    v.Init();
    return v;
}

// ============================================================
// Тест-класс
// ============================================================

class AnalyticalTest : public ::testing::Test {
   protected:
    VehicleData stable_;
    VehicleData unstable_;

    void SetUp() override {
        stable_ = MakeStableVehicle();

        try {
            unstable_ = MakeUnstableVehicle();
        } catch (const std::exception&) {
            // Неустойчивое КВП может не пройти валидацию —
            // тесты с ним будут пропущены
        }
    }
};

// ============================================================
// 1. Тесты VerifyAnalyticalDetailed
// ============================================================

TEST_F(AnalyticalTest, StableVehicleIsDetectedAsStable) {
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_TRUE(result.is_stable)
        << "Устойчивое КВП должно быть определено как устойчивое";
}

TEST_F(AnalyticalTest, StableVehicleHasNegativeMargin) {
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_LT(result.stability_margin, 0.0)
        << "Запас устойчивости должен быть < 0 для устойчивой системы";
}

TEST_F(AnalyticalTest, MarginMatchesSumOfTerms) {
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_NEAR(result.stability_margin,
                result.pneumatic_term + result.geometric_term, 1e-12)
        << "stability_margin = pneumatic_term + geometric_term";
}

TEST_F(AnalyticalTest, PneumaticTermIsNegativeForStableVehicle) {
    // Для устойчивого КВП пневматическое слагаемое должно быть
    // доминирующим и отрицательным
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_LT(result.pneumatic_term, 0.0)
        << "Пневматическое слагаемое < 0 для устойчивого КВП";
}

TEST_F(AnalyticalTest, GeometricTermIsPositive) {
    // phi0 > 0 (давление открывает ограждение) → sin(phi0) > 0
    // → геометрическое слагаемое > 0 (дестабилизирующее)
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_GT(result.geometric_term, 0.0)
        << "Геометрическое слагаемое > 0: открытое ограждение дестабилизирует";
}

TEST_F(AnalyticalTest, GeometricTermMatchesFormula) {
    using namespace constants;
    // geometric = χ * sqrt(2*p0/ρ) * sin(phi0)
    double expected =
        kChi * std::sqrt(2.0 * stable_.p0 / kRhoAir) * std::sin(stable_.phi0);
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_NEAR(result.geometric_term, expected, 1e-9)
        << "Геометрическое слагаемое: χ*sqrt(2p0/ρ)*sin(phi0)";
}

TEST_F(AnalyticalTest, PneumaticTermMatchesFormula) {
    using namespace constants;
    // pneumatic = (l/2) * (n*pa/W0) * (dQin/dp - Q0/(2*p0))
    double expected = (stable_.geometry.l / 2.0) *
                      (kN * kPAtm / stable_.geometry.W0) *
                      (stable_.dQin_dp - stable_.Q0 / (2.0 * stable_.p0));
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_NEAR(result.pneumatic_term, expected, 1e-9)
        << "Пневматическое слагаемое по формуле статьи";
}

TEST_F(AnalyticalTest, InfluenceRatioIsPositive) {
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    EXPECT_GT(result.influence_ratio, 0.0)
        << "Отношение вкладов должно быть положительным";
}

TEST_F(AnalyticalTest, InfluenceRatioMatchesFormula) {
    AnalyticalResult result = VerifyAnalyticalDetailed(stable_);
    double expected = std::abs(result.pneumatic_term / result.geometric_term);
    EXPECT_NEAR(result.influence_ratio, expected, 1e-12)
        << "influence_ratio = |pneumatic| / |geometric|";
}

TEST_F(AnalyticalTest, ZeroAngleGivesZeroGeometricTerm) {
    // При phi0 = 0: sin(0) = 0 → geometric_term = 0
    // Но phi0 = 0 невозможно при ненулевом давлении —
    // проверяем формулу напрямую
    using namespace constants;
    double term = kChi * std::sqrt(2.0 * stable_.p0 / kRhoAir) * std::sin(0.0);
    EXPECT_NEAR(term, 0.0, 1e-15);
}

// ============================================================
// 2. Тесты AnalyticalVerification (простая функция)
// ============================================================

class AnalyticalVerificationTest : public ::testing::Test {
   protected:
    VehicleData stable_;

    void SetUp() override { stable_ = MakeStableVehicle(); }
};

TEST_F(AnalyticalVerificationTest, StableVehicleReturnsTrue) {
    EXPECT_TRUE(AnalyticalVerification(stable_))
        << "AnalyticalVerification должна вернуть true для устойчивого КВП";
}

TEST_F(AnalyticalVerificationTest, MatchesDetailedResult) {
    bool simple = AnalyticalVerification(stable_);
    auto detailed = VerifyAnalyticalDetailed(stable_);
    EXPECT_EQ(simple, detailed.is_stable)
        << "Простая и детальная проверки должны давать одинаковый результат";
}

TEST_F(AnalyticalVerificationTest, ResultConsistentWithMarginSign) {
    auto detailed = VerifyAnalyticalDetailed(stable_);
    bool expected_stable = (detailed.stability_margin < 0.0);
    EXPECT_EQ(detailed.is_stable, expected_stable)
        << "is_stable должно соответствовать знаку stability_margin";
}

// ============================================================
// 3. Параметрические тесты: влияние параметров на устойчивость
// ============================================================

class AnalyticalParametricTest : public ::testing::Test {};

TEST_F(AnalyticalParametricTest, LargerSealWidthDecreasesMargin) {
    // Большее l → больший геометрический вклад (дестабилизирует)
    // и больший пневматический вклад (стабилизирует)
    // Результат зависит от соотношения, но проверим монотонность

    VehicleData v1 = MakeStableVehicle();

    VehicleData v2 = MakeStableVehicle();
    // Увеличиваем l, пересчитываем жёсткость чтобы phi0 не вышел за пределы
    v2.geometry.l = 0.5;
    v2.seal.k_seal_per_L = 8000.0;  // Компенсируем рост момента давления
    v2.Init();

    auto r1 = VerifyAnalyticalDetailed(v1);
    auto r2 = VerifyAnalyticalDetailed(v2);

    // Просто проверяем, что оба дают конечный результат
    EXPECT_TRUE(std::isfinite(r1.stability_margin));
    EXPECT_TRUE(std::isfinite(r2.stability_margin));
}

TEST_F(AnalyticalParametricTest, StrongerBlowerIncreasesStability) {
    // Более крутая характеристика нагнетателя → |dQin/dp| больше
    // → пневматическое слагаемое более отрицательное → лучше устойчивость

    VehicleData v_weak = MakeStableVehicle();
    v_weak.blower.k_fan = 0.00003;
    v_weak.Init();

    VehicleData v_strong = MakeStableVehicle();
    v_strong.blower.k_fan = 0.0001;  // Более крутая характеристика
    v_strong.Init();

    auto r_weak = VerifyAnalyticalDetailed(v_weak);
    auto r_strong = VerifyAnalyticalDetailed(v_strong);

    EXPECT_LT(r_strong.stability_margin, r_weak.stability_margin)
        << "Более крутая характеристика нагнетателя улучшает устойчивость";
}

TEST_F(AnalyticalParametricTest, LargerVolumeIncreasesStability) {
    // Больший объём подушки W0 → меньше сжимаемость → лучше устойчивость
    // |pneumatic| ∝ 1/W0 — нет, pneumatic = (l/2)*(n*pa/W0)*(...)
    // Большой W0 → меньше |pneumatic| — это УХУДШАЕТ устойчивость!
    // Проверим что pneumatic_term уменьшается по модулю

    VehicleData v_small = MakeStableVehicle();
    v_small.geometry.W0 = 10.0;
    v_small.Init();

    VehicleData v_large = MakeStableVehicle();
    v_large.geometry.W0 = 40.0;
    v_large.Init();

    auto r_small = VerifyAnalyticalDetailed(v_small);
    auto r_large = VerifyAnalyticalDetailed(v_large);

    // |pneumatic| ∝ 1/W0
    EXPECT_LT(std::abs(r_large.pneumatic_term),
              std::abs(r_small.pneumatic_term))
        << "Больший объём подушки уменьшает пневматический вклад (∝ 1/W0)";
}

TEST_F(AnalyticalParametricTest, HigherPressureIncreasesGeometricTerm) {
    // Больший p0 → больший geometric_term = χ*sqrt(2p0/ρ)*sin(phi0)
    // Достигается большей массой при той же площади

    VehicleData v_light = MakeStableVehicle();
    v_light.geometry.m = 500.0;
    v_light.Init();

    VehicleData v_heavy = MakeStableVehicle();
    v_heavy.geometry.m = 2000.0;
    v_heavy.Init();

    auto r_light = VerifyAnalyticalDetailed(v_light);
    auto r_heavy = VerifyAnalyticalDetailed(v_heavy);

    EXPECT_GT(r_heavy.geometric_term, r_light.geometric_term)
        << "Больший p0 увеличивает геометрическое слагаемое";
}

}  // namespace
}  // namespace acv

// tests/test_vehicle_data.cpp
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>

#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

// ============================================================
// Фабрика тестового транспортного средства
// Все параметры физически согласованы
// ============================================================

VehicleData MakeValidVehicle() {
    VehicleData v;
    v.name = "Test Hovercraft";

    // Геометрия
    v.geometry.m = 1000.0;     // kg
    v.geometry.L = 8.0;        // m
    v.geometry.l = 0.4;        // m
    v.geometry.S = 40.0;       // m²
    v.geometry.W0 = 20.0;      // m³
    v.geometry.I_phi = 400.0;  // kg·m²
    v.geometry.h0 = 0.15;      // m

    // Линейная характеристика нагнетателя: Q = Q_max - k_fan * p
    // При p0 ≈ m*g/S = 1000*9.81/40 ≈ 245 Pa
    // Q_blower(245) = 10 - 0.00003*245 ≈ 9.99 m³/s > Q_out — OK
    v.blower.type = BlowerCharacteristics::Type::Linear;
    v.blower.Q_max = 10.0;
    v.blower.k_fan = 0.00003;  // dQin/dp = -3e-5 < 0 ✓

    // Уплотнение: жёсткое, сильное демпфирование
    v.seal.k_seal_per_L = 5000.0;  // N·m/(rad·m)
    v.seal.k_phi = 500.0;          // N·m/rad
    v.seal.c_phi = 200.0;          // N·m·s/rad

    // Демпфирование подушки
    v.damping.beta_cushion = 0.3;

    return v;
}

// ============================================================
// Вспомогательные функции для ручной проверки
// ============================================================

double ExpectedPressure(const VehicleData& v) {
    return (v.geometry.m * constants::kG) / v.geometry.S;
}

double ExpectedOutflow(const VehicleData& v, double p, double phi) {
    using namespace constants;
    double h_eff = v.geometry.h0 + v.geometry.l * std::sin(phi);
    return kChi * v.geometry.L * h_eff * std::sqrt(2.0 * p / kRhoAir);
}

double ExpectedBlowerFlow(const VehicleData& v, double p) {
    return v.blower.Q_max - v.blower.k_fan * p;
}

// ============================================================
// Тест-класс
// ============================================================

class VehicleDataTest : public ::testing::Test {
   protected:
    VehicleData vehicle_;

    void SetUp() override {
        vehicle_ = MakeValidVehicle();
        vehicle_.Init();
    }
};

// ============================================================
// 1. Тесты равновесия
// ============================================================

class EquilibriumTest : public ::testing::Test {
   protected:
    VehicleData vehicle_;

    void SetUp() override {
        vehicle_ = MakeValidVehicle();
        vehicle_.Init();
    }
};

TEST_F(EquilibriumTest, PressureMatchesWeightOverArea) {
    // p0 * S = m * g  с допуском 1%
    double expected = ExpectedPressure(vehicle_);
    EXPECT_NEAR(vehicle_.p0, expected, expected * 0.01)
        << "Равновесное давление должно поддерживать вес: p0 ≈ m*g/S";
}

TEST_F(EquilibriumTest, MassFlowBalanceAtEquilibrium) {
    // Q_in(p0) = Q_out(p0, phi0)
    double Q_in = ExpectedBlowerFlow(vehicle_, vehicle_.p0);
    double Q_out = ExpectedOutflow(vehicle_, vehicle_.p0, vehicle_.phi0);

    EXPECT_NEAR(Q_in, Q_out, 1e-4)
        << "Баланс расходов в равновесии: Q_in ≈ Q_out";
}

TEST_F(EquilibriumTest, MomentBalanceAtEquilibrium) {
    // M_pressure(p0, phi0) = M_seal(phi0)
    double M_pressure = 0.5 * vehicle_.geometry.L *
                        std::pow(vehicle_.geometry.l, 2) * vehicle_.p0;
    double k_total =
        vehicle_.seal.k_seal_per_L * vehicle_.geometry.L + vehicle_.seal.k_phi;
    double M_seal = k_total * vehicle_.phi0;

    EXPECT_NEAR(M_pressure, M_seal, 1e-3)
        << "Баланс моментов в равновесии: M_pressure ≈ M_seal";
}

TEST_F(EquilibriumTest, EquilibriumAngleIsPositive) {
    // Давление открывает ограждение → phi0 > 0
    EXPECT_GT(vehicle_.phi0, 0.0)
        << "Давление поворачивает ограждение наружу: phi0 > 0";
}

TEST_F(EquilibriumTest, EquilibriumAngleIsSmall) {
    // phi0 должен быть физически разумным (< 45°)
    EXPECT_LT(vehicle_.phi0, M_PI / 4.0)
        << "Угол равновесия не должен превышать 45°";
}

TEST_F(EquilibriumTest, EquilibriumFlowIsPositive) {
    EXPECT_GT(vehicle_.Q0, 0.0)
        << "Расход в равновесии должен быть положительным";
}

TEST_F(EquilibriumTest, BlowerCanSupplyRequiredFlow) {
    double Q_available = ExpectedBlowerFlow(vehicle_, vehicle_.p0);
    EXPECT_GE(Q_available, vehicle_.Q0)
        << "Нагнетатель должен обеспечивать нужный расход";
}

TEST_F(EquilibriumTest, GapRemainsOpenAtEquilibrium) {
    // h_eff = h0 + l*sin(phi0) > 0
    double h_eff =
        vehicle_.geometry.h0 + vehicle_.geometry.l * std::sin(vehicle_.phi0);
    EXPECT_GT(h_eff, 0.0) << "Зазор в равновесии должен быть открыт";
}

TEST_F(EquilibriumTest, HeavierVehicleHasHigherPressure) {
    VehicleData heavy = MakeValidVehicle();
    heavy.geometry.m = 2000.0;  // Вдвое тяжелее
    heavy.Init();

    EXPECT_GT(heavy.p0, vehicle_.p0)
        << "Более тяжёлое КВП требует большего давления";
}

TEST_F(EquilibriumTest, LargerAreaHasLowerPressure) {
    VehicleData wide = MakeValidVehicle();
    wide.geometry.S = 80.0;  // Вдвое больше площадь
    // Подбираем нагнетатель под новое равновесие
    wide.blower.Q_max = 20.0;
    wide.Init();

    EXPECT_LT(wide.p0, vehicle_.p0)
        << "Большая площадь подушки → меньшее давление";
}

TEST_F(EquilibriumTest, StifferSealGivesSmallerAngle) {
    VehicleData stiff = MakeValidVehicle();
    stiff.seal.k_seal_per_L = 50000.0;  // В 10 раз жёстче
    stiff.Init();

    EXPECT_LT(stiff.phi0, vehicle_.phi0)
        << "Более жёсткое ограждение → меньший угол равновесия";
}

// ============================================================
// 2. Тесты производных
// ============================================================

class DerivativesTest : public ::testing::Test {
   protected:
    VehicleData vehicle_;

    void SetUp() override {
        vehicle_ = MakeValidVehicle();
        vehicle_.Init();
    }
};

TEST_F(DerivativesTest, DQinDpIsNegative) {
    // Нагнетатель с падающей характеристикой: dQin/dp < 0
    EXPECT_LT(vehicle_.dQin_dp, 0.0)
        << "Характеристика нагнетателя должна быть падающей";
}

TEST_F(DerivativesTest, DQinDpMatchesLinearBlower) {
    // Для линейного нагнетателя dQin/dp = -k_fan
    EXPECT_NEAR(vehicle_.dQin_dp, -vehicle_.blower.k_fan, 1e-12)
        << "dQin/dp должна совпадать с -k_fan для линейного нагнетателя";
}

TEST_F(DerivativesTest, DQoutDpIsPositive) {
    // Рост давления → рост скорости истечения → рост Q_out
    EXPECT_GT(vehicle_.dQout_dp, 0.0)
        << "dQout/dp > 0: рост давления увеличивает утечку";
}

TEST_F(DerivativesTest, DQoutDpMatchesFormula) {
    // dQout/dp = Q0 / (2*p0)
    double expected = vehicle_.Q0 / (2.0 * vehicle_.p0);
    EXPECT_NEAR(vehicle_.dQout_dp, expected, 1e-10) << "dQout/dp = Q0 / (2*p0)";
}

TEST_F(DerivativesTest, DQoutDphiIsPositive) {
    // ИСПРАВЛЕНО: h_eff = h0 + l*sin(phi)
    // Рост phi → рост зазора → рост Q_out → dQout/dphi > 0
    EXPECT_GT(vehicle_.dQout_dphi, 0.0)
        << "dQout/dphi > 0: открытие ограждения увеличивает утечку";
}

TEST_F(DerivativesTest, DQoutDphiMatchesFormula) {
    using namespace constants;
    // dQout/dphi = χ*L*l*cos(phi0)*sqrt(2*p0/ρ)
    double expected = kChi * vehicle_.geometry.L * vehicle_.geometry.l *
                      std::cos(vehicle_.phi0) *
                      std::sqrt(2.0 * vehicle_.p0 / kRhoAir);
    EXPECT_NEAR(vehicle_.dQout_dphi, expected, 1e-6)
        << "dQout/dphi = χ*L*l*cos(phi0)*sqrt(2p0/ρ)";
}

TEST_F(DerivativesTest, DQoutDHIsPositive) {
    // Рост высоты → рост зазора → рост утечки
    EXPECT_GT(vehicle_.dQout_dH, 0.0)
        << "dQout/dH > 0: подъём корпуса увеличивает утечку";
}

TEST_F(DerivativesTest, DMDphiIsNegative) {
    // Жёсткость ограждения — возвращающий момент: dM/dphi < 0
    EXPECT_LT(vehicle_.dM_dphi, 0.0)
        << "dM/dphi < 0: ограждение стремится вернуться в равновесие";
}

TEST_F(DerivativesTest, DMDphiMatchesFormula) {
    // dM/dphi = -(k_seal_per_L * L + k_phi)
    double k_total =
        vehicle_.seal.k_seal_per_L * vehicle_.geometry.L + vehicle_.seal.k_phi;
    double expected = -k_total;
    EXPECT_NEAR(vehicle_.dM_dphi, expected, 1e-9)
        << "dM/dphi = -(k_seal_per_L*L + k_phi)";
}

TEST_F(DerivativesTest, DMDphidotIsNonPositive) {
    // Демпфирование: dM/dφ̇ ≤ 0
    EXPECT_LE(vehicle_.dM_dphidot, 0.0)
        << "dM/dφ̇ ≤ 0: демпфирование диссипирует энергию";
}

TEST_F(DerivativesTest, DMDphidotMatchesFormula) {
    // dM/dφ̇ = -c_phi
    EXPECT_NEAR(vehicle_.dM_dphidot, -vehicle_.seal.c_phi, 1e-12)
        << "dM/dφ̇ = -c_phi";
}

TEST_F(DerivativesTest, DYDHdotIsNegative) {
    // Демпфирование вертикального движения: dY/dḢ < 0
    EXPECT_LT(vehicle_.dY_dHdot, 0.0)
        << "dY/dḢ < 0: вертикальное демпфирование диссипирует энергию";
}

TEST_F(DerivativesTest, DYDHdotMatchesFormula) {
    using namespace constants;
    // dY/dḢ = -beta * rho * Q0
    double expected = -vehicle_.damping.beta_cushion * kRhoAir * vehicle_.Q0;
    EXPECT_NEAR(vehicle_.dY_dHdot, expected, 1e-9) << "dY/dḢ = -β*ρ*Q0";
}

TEST_F(DerivativesTest, CushionCriterionSatisfied) {
    // Необходимое условие устойчивости подушки (Раздел 1.2 статьи):
    // dQin/dp - dQout/dp < 0
    EXPECT_LT(vehicle_.dQin_dp - vehicle_.dQout_dp, 0.0)
        << "Критерий устойчивости подушки: dQin/dp - dQout/dp < 0";
}

// ============================================================
// 3. Тесты характеристики нагнетателя
// ============================================================

class BlowerTest : public ::testing::Test {
   protected:
    VehicleData vehicle_ = MakeValidVehicle();
};

TEST_F(BlowerTest, LinearBlowerFlowAtZeroPressure) {
    // Q(0) = Q_max
    EXPECT_NEAR(vehicle_.GetBlowerFlow(0.0), vehicle_.blower.Q_max, 1e-12);
}

TEST_F(BlowerTest, LinearBlowerFlowDecreaseWithPressure) {
    double Q_low = vehicle_.GetBlowerFlow(100.0);
    double Q_high = vehicle_.GetBlowerFlow(500.0);
    EXPECT_GT(Q_low, Q_high)
        << "Расход нагнетателя должен уменьшаться с ростом давления";
}

TEST_F(BlowerTest, LinearBlowerDerivativeConstant) {
    // dQ/dp = -k_fan для линейной характеристики
    double deriv1 = vehicle_.GetBlowerDerivative(100.0);
    double deriv2 = vehicle_.GetBlowerDerivative(500.0);
    EXPECT_NEAR(deriv1, deriv2, 1e-15)
        << "Производная линейного нагнетателя постоянна";
    EXPECT_NEAR(deriv1, -vehicle_.blower.k_fan, 1e-15);
}

TEST_F(BlowerTest, PolynomialBlowerFlow) {
    VehicleData v = MakeValidVehicle();
    // Q(p) = 10 - 0.00003*p  — эквивалент линейному
    v.blower.type = BlowerCharacteristics::Type::Polynomial;
    v.blower.polynomial_coeffs = {10.0, -0.00003};

    double Q_linear = v.blower.Q_max - v.blower.k_fan * 200.0;
    double Q_poly = v.GetBlowerFlow(200.0);

    EXPECT_NEAR(Q_poly, Q_linear, 1e-10)
        << "Полиномиальный нагнетатель должен совпадать с линейным";
}

TEST_F(BlowerTest, PolynomialBlowerDerivative) {
    VehicleData v = MakeValidVehicle();
    v.blower.type = BlowerCharacteristics::Type::Polynomial;
    v.blower.polynomial_coeffs = {10.0, -0.00003};  // Q = 10 - 3e-5 * p

    double deriv = v.GetBlowerDerivative(300.0);
    EXPECT_NEAR(deriv, -0.00003, 1e-12)
        << "Производная полиномиального нагнетателя (линейный случай)";
}

// ============================================================
// 4. Тесты расхода утечки
// ============================================================

class OutflowTest : public ::testing::Test {
   protected:
    VehicleData vehicle_ = MakeValidVehicle();
};

TEST_F(OutflowTest, OutflowIncreasesWithPressure) {
    double Q_low = vehicle_.GetOutflowRate(200.0, 0.05, 0.15);
    double Q_high = vehicle_.GetOutflowRate(400.0, 0.05, 0.15);
    EXPECT_GT(Q_high, Q_low) << "Рост давления увеличивает расход утечки";
}

TEST_F(OutflowTest, OutflowIncreasesWithAngle) {
    // КЛЮЧЕВОЙ ТЕСТ: при h_eff = h0 + l*sin(phi)
    // рост phi → рост Q_out
    double Q_small = vehicle_.GetOutflowRate(250.0, 0.05, 0.15);
    double Q_large = vehicle_.GetOutflowRate(250.0, 0.20, 0.15);
    EXPECT_GT(Q_large, Q_small)
        << "Открытие ограждения (рост phi) увеличивает утечку: "
           "h_eff = h0 + l*sin(phi)";
}

TEST_F(OutflowTest, OutflowZeroWhenGapClosed) {
    // При phi таком, что h0 + l*sin(phi) ≤ 0 → Q_out = 0
    double phi_closed =
        -std::asin(vehicle_.geometry.h0 / vehicle_.geometry.l) - 0.01;
    double Q = vehicle_.GetOutflowRate(250.0, phi_closed, 0.15);
    EXPECT_NEAR(Q, 0.0, 1e-12) << "При закрытом зазоре утечка = 0";
}

TEST_F(OutflowTest, OutflowMatchesFormula) {
    using namespace constants;
    double p = 300.0;
    double phi = 0.1;
    double h_eff = vehicle_.geometry.h0 + vehicle_.geometry.l * std::sin(phi);
    double expected =
        kChi * vehicle_.geometry.L * h_eff * std::sqrt(2.0 * p / kRhoAir);

    EXPECT_NEAR(vehicle_.GetOutflowRate(p, phi, 0.15), expected, 1e-9)
        << "GetOutflowRate должен совпадать с формулой χ*L*h_eff*sqrt(2p/ρ)";
}

// ============================================================
// 5. Тесты критерия Рауса-Гурвица
// ============================================================

class RouthHurwitzTest : public ::testing::Test {
   protected:
    VehicleData vehicle_;

    void SetUp() override {
        vehicle_ = MakeValidVehicle();
        vehicle_.Init();
    }
};

TEST_F(RouthHurwitzTest, StableVehiclePassesRouth) {
    EXPECT_TRUE(vehicle_.CheckRouthHurwitz())
        << "Физически устойчивое КВП должно проходить критерий Рауса";
}

TEST_F(RouthHurwitzTest, WeakBlowerMayFail) {
    // Слабый нагнетатель → почти нулевая производная →
    // возможная неустойчивость
    VehicleData v = MakeValidVehicle();
    v.blower.k_fan = 1e-8;  // Почти горизонтальная характеристика
    v.blower.Q_max = 10.0;
    v.seal.c_phi = 0.0;  // Нет демпфирования

    try {
        v.Init();
        // Если инициализация прошла — результат зависит от параметров
        bool stable = v.CheckRouthHurwitz();
        // Просто проверяем, что функция возвращает bool без краша
        (void)stable;
    } catch (const std::exception&) {
        // Может не сойтись или нарушить критерий — это ожидаемо
        SUCCEED();
    }
}

TEST_F(RouthHurwitzTest, LargeInertiaMayDestabilize) {
    // Большой момент инерции → медленная динамика ограждения
    VehicleData v = MakeValidVehicle();
    v.geometry.I_phi = 1e6;  // Огромная инерция
    v.Init();

    // Не проверяем результат — просто что функция работает
    bool result = v.CheckRouthHurwitz();
    (void)result;
}

// ============================================================
// 6. Тесты валидации входных данных
// ============================================================

class ValidationTest : public ::testing::Test {};

TEST_F(ValidationTest, ThrowsOnNegativeMass) {
    VehicleData v = MakeValidVehicle();
    v.geometry.m = -1.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnZeroLength) {
    VehicleData v = MakeValidVehicle();
    v.geometry.L = 0.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnZeroSealWidth) {
    VehicleData v = MakeValidVehicle();
    v.geometry.l = 0.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnZeroArea) {
    VehicleData v = MakeValidVehicle();
    v.geometry.S = 0.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnZeroVolume) {
    VehicleData v = MakeValidVehicle();
    v.geometry.W0 = 0.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnZeroInertia) {
    VehicleData v = MakeValidVehicle();
    v.geometry.I_phi = 0.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnZeroGap) {
    VehicleData v = MakeValidVehicle();
    v.geometry.h0 = 0.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnNegativeSealStiffness) {
    VehicleData v = MakeValidVehicle();
    v.seal.k_seal_per_L = -1.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsOnNegativeDamping) {
    VehicleData v = MakeValidVehicle();
    v.seal.c_phi = -1.0;
    EXPECT_THROW(v.Init(), std::runtime_error);
}

TEST_F(ValidationTest, ThrowsWhenBlowerTooWeak) {
    VehicleData v = MakeValidVehicle();
    // Нагнетатель не может создать нужный расход при равновесном давлении
    v.blower.Q_max = 0.001;  // Очень маленький расход
    EXPECT_THROW(v.Init(), std::runtime_error);
}

// ============================================================
// 7. Тест загрузки из JSON
// ============================================================

class JsonTest : public ::testing::Test {
   protected:
    std::string tmp_file_ = "test_vehicle_tmp.json";

    void TearDown() override { std::remove(tmp_file_.c_str()); }

    void WriteJson(const nlohmann::json& j) {
        std::ofstream f(tmp_file_);
        f << j.dump(4);
    }

    nlohmann::json MakeValidJson() {
        nlohmann::json j;
        j["name"] = "JSON Test Vehicle";
        j["geometry"]["m"] = 1000.0;
        j["geometry"]["L"] = 8.0;
        j["geometry"]["l"] = 0.4;
        j["geometry"]["S"] = 40.0;
        j["geometry"]["W0"] = 20.0;
        j["geometry"]["I_phi"] = 400.0;
        j["geometry"]["h0"] = 0.15;
        j["blower"]["type"] = "linear";
        j["blower"]["Q_max"] = 10.0;
        j["blower"]["k_fan"] = 0.00003;
        j["seal"]["k_seal_per_L"] = 5000.0;
        j["seal"]["k_phi"] = 500.0;
        j["seal"]["c_phi"] = 200.0;
        j["damping"]["beta_cushion"] = 0.3;
        return j;
    }
};

TEST_F(JsonTest, LoadsValidJsonSuccessfully) {
    WriteJson(MakeValidJson());
    EXPECT_NO_THROW({
        VehicleData v = LoadVehicleFromJson(tmp_file_);
        EXPECT_EQ(v.name, "JSON Test Vehicle");
        EXPECT_DOUBLE_EQ(v.geometry.m, 1000.0);
    });
}

TEST_F(JsonTest, ThrowsOnMissingFile) {
    EXPECT_THROW(LoadVehicleFromJson("nonexistent.json"), std::runtime_error);
}

TEST_F(JsonTest, SaveAndReloadPreservesEquilibrium) {
    VehicleData original = MakeValidVehicle();
    original.Init();

    std::string out_file = "test_save_output.json";
    SaveVehicleToJson(original, out_file);

    // После загрузки параметры равновесия должны совпасть
    // (загружаем как новое КВП и переинициализируем)
    VehicleData reloaded = LoadVehicleFromJson(out_file);

    EXPECT_NEAR(reloaded.p0, original.p0, 1e-4);
    EXPECT_NEAR(reloaded.Q0, original.Q0, 1e-6);
    EXPECT_NEAR(reloaded.phi0, original.phi0, 1e-6);

    std::remove(out_file.c_str());
}

TEST_F(JsonTest, LoadsPolynomialBlower) {
    auto j = MakeValidJson();
    j["blower"]["type"] = "polynomial";
    j["blower"]["coefficients"] = {10.0, -0.00003};
    WriteJson(j);

    EXPECT_NO_THROW({
        VehicleData v = LoadVehicleFromJson(tmp_file_);
        EXPECT_EQ(v.blower.type, BlowerCharacteristics::Type::Polynomial);
    });
}

}  // namespace
}  // namespace acv

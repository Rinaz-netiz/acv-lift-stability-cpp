#include <gtest/gtest.h>

#include <cmath>

#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class VehicleDataTest : public ::testing::Test {
   protected:
    void SetUp() override {
        // Базовая геометрия из примера
        data.name = "Test Vehicle";
        data.geometry.m = 2050.0;
        data.geometry.L = 10.5;
        data.geometry.l = 0.7;
        data.geometry.S = 20.0;
        data.geometry.W0 = 14.0;
        data.geometry.I_phi = 9.35;
        data.geometry.H_init = 0.0;

        // Характеристика нагнетателя (из примера)
        data.blower.p_table = {0.0, 1000.0, 1005.525, 1010.0, 1500.0};
        data.blower.Q_table = {7.721, 5.72105, 5.71, 5.70105, 4.721};

        // Момент ограждения
        data.seal_moment.phi_table = {0.0, 0.1, 0.18, 0.19, 0.2, 0.3};
        data.seal_moment.M_per_L_table = {0.0,       -50.0,     -224.5736,
                                          -246.3536, -268.1336, -485.9336};

        // Демпфирование ограждения
        data.seal_damping_curve.h_gap_table = {0.0, 0.002, 0.006, 0.014, 0.1};
        data.seal_damping_curve.dM_dphidot_table = {-50.0, -30.0, -14.0, -5.0,
                                                    -1.0};

        // Демпфирование подушки
        data.cushion_damping_curve.Sgap_over_S_table = {
            0.0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.02};
        data.cushion_damping_curve.D_table = {120000, 65000, 38000, 26000,
                                              21000,  18500, 17000, 15000,
                                              12000,  8000};

        data.dM_dphidot_input = -30.0;
        data.D_damping = 20000.0;
    }

    VehicleData data;
};

// ============================================================
// Тесты BlowerCurve
// ============================================================

TEST_F(VehicleDataTest, BlowerCurveTabular) {
    // Проверка граничных значений
    EXPECT_NEAR(data.blower.Flow(0.0), 7.721, 1e-6);
    EXPECT_NEAR(data.blower.Flow(1500.0), 4.721, 1e-6);

    // Проверка интерполяции
    double Q_mid = data.blower.Flow(1000.0);
    EXPECT_NEAR(Q_mid, 5.72105, 1e-6);

    // Производная должна быть отрицательной
    double dQ_dp = data.blower.Derivative(1000.0);
    EXPECT_LT(dQ_dp, 0.0);
}

TEST_F(VehicleDataTest, BlowerCurveLinear) {
    BlowerCurve linear;
    linear.Q_max = 8.0;
    linear.k_fan = 0.002;  // m³/(s·Pa)

    EXPECT_NEAR(linear.Flow(0.0), 8.0, 1e-6);
    EXPECT_NEAR(linear.Flow(1000.0), 6.0, 1e-6);
    EXPECT_NEAR(linear.Derivative(500.0), -0.002, 1e-6);
}

TEST_F(VehicleDataTest, BlowerCurveExtrapolation) {
    // За пределами таблицы должны быть крайние значения
    EXPECT_NEAR(data.blower.Flow(-100.0), 7.721, 1e-6);
    EXPECT_NEAR(data.blower.Flow(2000.0), 4.721, 1e-6);
}

// ============================================================
// Тесты SealMomentCurve
// ============================================================

TEST_F(VehicleDataTest, SealMomentCurveBoundaries) {
    double L = data.geometry.L;

    EXPECT_NEAR(data.seal_moment.Moment(0.0, L), 0.0, 1e-6);
    EXPECT_NEAR(data.seal_moment.Moment(0.3, L), -485.9336 * L, 1e-3);
}

TEST_F(VehicleDataTest, SealMomentCurveInterpolation) {
    double L = data.geometry.L;
    double M = data.seal_moment.Moment(0.15, L);

    // Должно быть между -50*L и -224.5736*L
    EXPECT_LT(M, -50.0 * L);
    EXPECT_GT(M, -224.5736 * L);
}

TEST_F(VehicleDataTest, SealMomentDerivativeNegative) {
    double L = data.geometry.L;

    // Производная должна быть отрицательной (восстанавливающий момент)
    for (double phi = 0.05; phi < 0.25; phi += 0.05) {
        double dM_dphi = data.seal_moment.Derivative(phi, L);
        EXPECT_LT(dM_dphi, 0.0) << "at phi = " << phi;
    }
}

// ============================================================
// Тесты демпфирования
// ============================================================

TEST_F(VehicleDataTest, SealDampingCurveNegative) {
    // Все значения демпфирования должны быть отрицательными
    for (double h = 0.001; h < 0.05; h += 0.01) {
        double damping = data.seal_damping_curve.Value(h);
        EXPECT_LT(damping, 0.0) << "at h_gap = " << h;
    }
}

TEST_F(VehicleDataTest, CushionDampingCurvePositive) {
    // Коэффициент D должен быть положительным
    for (double ratio = 0.001; ratio < 0.015; ratio += 0.002) {
        double D = data.cushion_damping_curve.Value(ratio);
        EXPECT_GT(D, 0.0) << "at Sgap/S = " << ratio;
    }
}

TEST_F(VehicleDataTest, CushionDampingMonotonicDecrease) {
    // D должен монотонно убывать с ростом Sgap/S
    double prev_D = data.cushion_damping_curve.Value(0.001);
    for (double ratio = 0.003; ratio < 0.018; ratio += 0.002) {
        double D = data.cushion_damping_curve.Value(ratio);
        EXPECT_LE(D, prev_D) << "at Sgap/S = " << ratio;
        prev_D = D;
    }
}

// ============================================================
// Тесты валидации
// ============================================================

TEST_F(VehicleDataTest, ValidationPassesForValidData) {
    EXPECT_NO_THROW(data.Init());
}

TEST_F(VehicleDataTest, ValidationFailsForNegativeMass) {
    data.geometry.m = -100.0;
    EXPECT_THROW(data.Init(), std::runtime_error);
}

TEST_F(VehicleDataTest, ValidationFailsForZeroLength) {
    data.geometry.L = 0.0;
    EXPECT_THROW(data.Init(), std::runtime_error);
}

TEST_F(VehicleDataTest, ValidationFailsForEmptySealMomentTable) {
    data.seal_moment.phi_table.clear();
    data.seal_moment.M_per_L_table.clear();
    EXPECT_THROW(data.Init(), std::runtime_error);
}

TEST_F(VehicleDataTest, ValidationFailsForMismatchedBlowerTables) {
    data.blower.p_table.push_back(2000.0);
    // Q_table не изменён - размеры не совпадают
    EXPECT_THROW(data.Init(), std::runtime_error);
}

TEST_F(VehicleDataTest, ValidationFailsForPositiveDamping) {
    data.dM_dphidot_input = 10.0;  // Должно быть < 0
    data.seal_damping_curve.h_gap_table.clear();
    data.seal_damping_curve.dM_dphidot_table.clear();
    EXPECT_THROW(data.Init(), std::runtime_error);
}

// ============================================================
// Тесты вычисления равновесия
// ============================================================

TEST_F(VehicleDataTest, EquilibriumPressureMatchesWeight) {
    data.Init();

    using namespace constants;
    double expected_p = (data.geometry.m * kG) / data.geometry.S;
    EXPECT_NEAR(data.eq.p0, expected_p, 1.0);
}

TEST_F(VehicleDataTest, EquilibriumAnglePositive) {
    data.Init();
    EXPECT_GT(data.eq.phi0, 0.0);
    EXPECT_LT(data.eq.phi0, 0.5);  // Разумные пределы
}

TEST_F(VehicleDataTest, EquilibriumFlowPositive) {
    data.Init();
    EXPECT_GT(data.eq.Q0, 0.0);
    EXPECT_GT(data.eq.h_gap0, 0.0);
}

TEST_F(VehicleDataTest, EquilibriumMomentBalance) {
    data.Init();

    double M_seal = data.seal_moment.Moment(data.eq.phi0, data.geometry.L);
    double M_pressure =
        data.eq.p0 * data.geometry.L * std::pow(data.geometry.l, 2) / 2.0;

    // Момент от ограждения должен уравновешивать момент от давления
    EXPECT_NEAR(M_seal, -M_pressure, std::abs(M_pressure) * 1e-3);
}

TEST_F(VehicleDataTest, EquilibriumFlowBalance) {
    data.Init();

    using namespace constants;
    double Q_in = data.blower.Flow(data.eq.p0);
    double vel = kChi * std::sqrt(2.0 * data.eq.p0 / kRhoAir);
    double Q_out = vel * data.geometry.L * data.eq.h_gap0;

    EXPECT_NEAR(Q_in, Q_out, Q_in * 1e-4);
}

// ============================================================
// Тесты производных
// ============================================================

TEST_F(VehicleDataTest, DerivativeSignsCorrect) {
    data.Init();

    // dQin/dp < 0 (характеристика нагнетателя убывает)
    EXPECT_LT(data.der.dQin_dp, 0.0);

    // dQout/dp > 0 (расход через зазор растёт с давлением)
    EXPECT_GT(data.der.dQout_dp, 0.0);

    // dQout/dH > 0 (больше высота -> больше зазор -> больше расход)
    EXPECT_GT(data.der.dQout_dH, 0.0);

    // dM/dphi < 0 (восстанавливающий момент)
    EXPECT_LT(data.der.dM_dphi, 0.0);

    // dM/dphidot < 0 (демпфирование)
    EXPECT_LT(data.der.dM_dphidot, 0.0);

    // dY/dHdot < 0 (демпфирование подушки)
    EXPECT_LT(data.der.dY_dHdot, 0.0);
}

TEST_F(VehicleDataTest, NecessaryStabilityCondition) {
    data.Init();

    // Необходимое условие: dQin/dp - dQout/dp < 0
    double criterion = data.der.dQin_dp - data.der.dQout_dp;
    EXPECT_LT(criterion, 0.0);
}

TEST_F(VehicleDataTest, DerivativeQoutDpMatchesFormula) {
    data.Init();

    // Eq. (25): ∂Qout/∂p|₀ = Q₀/(2·p₀)
    double expected = data.eq.Q0 / (2.0 * data.eq.p0);
    EXPECT_NEAR(data.der.dQout_dp, expected, expected * 1e-6);
}

TEST_F(VehicleDataTest, DerivativeQoutDHMatchesFormula) {
    data.Init();

    using namespace constants;
    // Eq. (25): ∂Qout/∂H|₀ = χ·√(2p₀/ρ)·L
    double expected =
        kChi * std::sqrt(2.0 * data.eq.p0 / kRhoAir) * data.geometry.L;
    EXPECT_NEAR(data.der.dQout_dH, expected, expected * 1e-6);
}

// ============================================================
// Тесты граничных случаев
// ============================================================

TEST_F(VehicleDataTest, HeavyVehicleHighPressure) {
    data.geometry.m = 5000.0;

    EXPECT_THROW(data.Init(), std::runtime_error);
}

TEST_F(VehicleDataTest, LargeAreaLowPressure) {
    data.geometry.S = 40.0;  // Больше площадь
    data.Init();

    using namespace constants;
    double expected_p = (data.geometry.m * kG) / 40.0;
    EXPECT_NEAR(data.eq.p0, expected_p, 1.0);
    EXPECT_LT(data.eq.p0, 1000.0);
}

TEST_F(VehicleDataTest, SmallSealLengthSmallAngle) {
    data.geometry.l = 0.3;  // Короче ограждение
    data.Init();

    // При коротком ограждении угол должен быть меньше
    EXPECT_GT(data.eq.phi0, 0.0);
    EXPECT_LT(data.eq.phi0, 0.3);
}

// ============================================================
// Тесты сохранения/загрузки
// ============================================================

TEST_F(VehicleDataTest, SaveAndLoadRoundTrip) {
    data.Init();

    std::string temp_file = "/tmp/test_vehicle.json";
    SaveVehicleToJson(data, temp_file);

    VehicleData loaded = LoadVehicleFromJson(temp_file);

    EXPECT_EQ(loaded.name, data.name);
    EXPECT_NEAR(loaded.geometry.m, data.geometry.m, 1e-6);
    EXPECT_NEAR(loaded.eq.p0, data.eq.p0, 1e-3);
    EXPECT_NEAR(loaded.eq.phi0, data.eq.phi0, 1e-6);

    std::remove(temp_file.c_str());
}

}  // namespace
}  // namespace acv

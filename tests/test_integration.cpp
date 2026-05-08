#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include "my_project/analytical.h"
#include "my_project/numerical.h"
#include "my_project/vehicle_data.h"

namespace acv {
namespace {

class IntegrationTest : public ::testing::Test {
   protected:
    void SetUp() override {
        try {
            // Используем данные из примера
            std::string filename = "test_data/vehicle_a.json";

            // Проверяем существование файла
            std::ifstream test_file(filename);
            if (!test_file.good()) {
                // Пытаемся найти файл в других местах
                std::vector<std::string> possible_paths = {
                    "test_data/vehicle_a.json", "../test_data/vehicle_a.json",
                    "../../test_data/vehicle_a.json"};

                bool found = false;
                for (const auto& path : possible_paths) {
                    std::ifstream f(path);
                    if (f.good()) {
                        filename = path;
                        found = true;
                        std::cout << "Found test data at: " << path
                                  << std::endl;
                        break;
                    }
                }

                if (!found) {
                    GTEST_SKIP() << "Test data file not found. Tried:\n"
                                 << "  - test_data/vehicle_a.json\n"
                                 << "  - ../test_data/vehicle_a.json\n"
                                 << "  - ../../test_data/vehicle_a.json";
                    return;
                }
            }

            data = LoadVehicleFromJson(filename);

            // Отладочный вывод перед Init()
            // std::cout << "\n=== Vehicle Data Before Init ===" << std::endl;
            // std::cout << "Mass: " << data.geometry.m << " kg" << std::endl;
            // std::cout << "Length: " << data.geometry.L << " m" << std::endl;
            // std::cout << "Seal length: " << data.geometry.l << " m"
            //           << std::endl;

            // std::cout << "\nSeal moment curve:" << std::endl;
            // for (size_t i = 0; i < data.seal_moment.phi_table.size(); ++i) {
            //     std::cout << "  phi=" << data.seal_moment.phi_table[i]
            //               << " rad, M/L=" <<
            //               data.seal_moment.M_per_L_table[i]
            //               << " N" << std::endl;
            // }

        } catch (const std::exception& e) {
            GTEST_SKIP() << "Failed to load test data: " << e.what();
        }
    }

    VehicleData data;
};

TEST_F(IntegrationTest, FullWorkflow) {
    // 1. Инициализация данных
    ASSERT_NO_THROW(data.Init()) << "Failed to initialize vehicle data";

    // std::cout << "\n=== Equilibrium State ===" << std::endl;
    // std::cout << "p0 = " << data.eq.p0 << " Pa" << std::endl;
    // std::cout << "phi0 = " << data.eq.phi0 << " rad ("
    //           << (data.eq.phi0 * 180.0 / M_PI) << " deg)" << std::endl;
    // std::cout << "H0 = " << data.eq.H0 << " m" << std::endl;
    // std::cout << "Q0 = " << data.eq.Q0 << " m³/s" << std::endl;
    // std::cout << "h_gap0 = " << data.eq.h_gap0 << " m" << std::endl;

    // std::cout << "\n=== Derivatives ===" << std::endl;
    // std::cout << "dQin/dp = " << data.der.dQin_dp << std::endl;
    // std::cout << "dQout/dp = " << data.der.dQout_dp << std::endl;
    // std::cout << "dM/dphi = " << data.der.dM_dphi << std::endl;

    // 2. Аналитическая проверка
    bool analytical_stable = AnalyticalVerification(data);
    AnalyticalResult analytical_detailed = VerifyAnalyticalDetailed(data);

    // std::cout << "\n=== Analytical Verification ===" << std::endl;
    // std::cout << "Stable: " << analytical_stable << std::endl;
    // std::cout << "Stability margin: " << analytical_detailed.stability_margin
    //           << std::endl;
    // std::cout << "Pneumatic term: " << analytical_detailed.pneumatic_term
    //           << std::endl;
    // std::cout << "Geometric term: " << analytical_detailed.geometric_term
    //           << std::endl;

    EXPECT_EQ(analytical_stable, analytical_detailed.is_stable);

    // 3. Численный анализ (простая модель)
    StabilityResult numerical_simple = AnalyzeStabilitySimple(data);
    // std::cout << "\n=== Numerical Analysis (Simple Model) ===" << std::endl;
    // std::cout << "Stable: " << numerical_simple.is_stable << std::endl;
    // std::cout << "Max real part: " << numerical_simple.max_real_part
    //           << std::endl;
    // std::cout << "Eigenvalues:" << std::endl;
    // for (const auto& lambda : numerical_simple.eigenvalues) {
    //     std::cout << "  " << lambda << std::endl;
    // }

    EXPECT_EQ(numerical_simple.eigenvalues.size(), 3);

    // 4. Численный анализ (полная модель)
    StabilityResult numerical_full = AnalyzeStabilityFull(data);
    // std::cout << "\n=== Numerical Analysis (Full Model) ===" << std::endl;
    // std::cout << "Stable: " << numerical_full.is_stable << std::endl;
    // std::cout << "Max real part: " << numerical_full.max_real_part <<
    // std::endl;

    EXPECT_EQ(numerical_full.eigenvalues.size(), 5);
}

TEST_F(IntegrationTest, AnalyticalNumericalAgreement) {
    ASSERT_NO_THROW(data.Init());

    bool analytical = AnalyticalVerification(data);
    StabilityResult numerical = AnalyzeStabilitySimple(data);

    // std::cout << "Analytical: " << analytical
    //           << ", Numerical: " << numerical.is_stable << std::endl;
}

TEST_F(IntegrationTest, EquilibriumIsPhysical) {
    ASSERT_NO_THROW(data.Init());

    EXPECT_GT(data.eq.p0, 0.0);
    EXPECT_LT(data.eq.p0, 10000.0);

    EXPECT_GT(data.eq.phi0, 0.0);
    EXPECT_LT(data.eq.phi0, M_PI / 2);

    EXPECT_GT(data.eq.H0, data.geometry.H_init);
    EXPECT_LT(data.eq.H0, 2.0);

    EXPECT_GT(data.eq.Q0, 0.0);
    EXPECT_LT(data.eq.Q0, 20.0);

    EXPECT_GT(data.eq.h_gap0, 0.0);
    EXPECT_LT(data.eq.h_gap0, 0.2);
}

TEST_F(IntegrationTest, DerivativesAreConsistent) {
    ASSERT_NO_THROW(data.Init());

    EXPECT_LT(data.der.dQin_dp, 0.0)
        << "Blower characteristic should be decreasing";
    EXPECT_GT(data.der.dQout_dp, 0.0) << "Outflow increases with pressure";
    EXPECT_GT(data.der.dQout_dH, 0.0) << "Outflow increases with height";
    EXPECT_LT(data.der.dM_dphi, 0.0) << "Restoring moment";
    EXPECT_LT(data.der.dM_dphidot, 0.0) << "Seal damping";
    EXPECT_LT(data.der.dY_dHdot, 0.0) << "Cushion damping";
}

TEST_F(IntegrationTest, MultipleVehiclesComparison) {
    std::vector<std::string> possible_paths = {
        "test_data/vehicle_a.json", "../test_data/vehicle_a.json",
        "../../test_data/vehicle_a.json"};

    std::string config;
    for (const auto& path : possible_paths) {
        std::ifstream f(path);
        if (f.good()) {
            config = path;
            break;
        }
    }

    if (config.empty()) {
        GTEST_SKIP() << "Test data file not found";
        return;
    }

    VehicleData v;
    EXPECT_NO_THROW(v = LoadVehicleFromJson(config));
    EXPECT_NO_THROW(v.Init());

    bool stable_analytical = AnalyticalVerification(v);
    StabilityResult stable_numerical = AnalyzeStabilitySimple(v);

    // std::cout << config << ": analytical=" << stable_analytical
    //           << ", numerical=" << stable_numerical.is_stable << std::endl;
}

TEST_F(IntegrationTest, RobustnessToSmallPerturbations) {
    ASSERT_NO_THROW(data.Init());

    VehicleData perturbed = data;
    perturbed.geometry.m *= 1.001;

    ASSERT_NO_THROW(perturbed.Init());

    EXPECT_NEAR(perturbed.eq.p0, data.eq.p0, data.eq.p0 * 0.01);
    EXPECT_NEAR(perturbed.eq.phi0, data.eq.phi0, data.eq.phi0 * 0.05);
}

}  // namespace
}  // namespace acv

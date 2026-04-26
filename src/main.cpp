// src/main.cpp
#include <iostream>

#include "my_project/analytical.h"
#include "my_project/chec.h"
#include "my_project/numerical.h"

int main(int argc, char* argv[]) {
    auto v = acv::loadVehicleFromJson(argv[1]);

    // Аналитическая модель
    bool stable_analytical = acv::analyticalVerification(v);

    // Численная модель
    auto coeffs = acv::computeFullCharacteristicPolynomial(v);
    bool stable_numerical = acv::isStableNumerical(coeffs);

    std::cout << "Analytical: " << (stable_analytical ? "STABLE" : "UNSTABLE")
              << "\n";
    std::cout << "Numerical:  " << (stable_numerical ? "STABLE" : "UNSTABLE")
              << "\n";

    acv::StabilitySolver::Analyze(v);

    try {
        auto res = acv::AnalyzeStability(v);

        std::cout << "--- Анализ устойчивости ---" << std::endl;
        std::cout << "Статус: " << (res.is_stable ? "УСТОЙЧИВА" : "НЕУСТОЙЧИВА")
                  << std::endl;
        std::cout << "Макс. вещественная часть: " << res.max_real_part
                  << std::endl;

        for (const auto& ev : res.eigenvalues) {
            std::cout << "Корень: " << ev.real() << " + " << ev.imag() << "i"
                      << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
    }
}

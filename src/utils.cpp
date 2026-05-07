// src/utils.cpp
#include <iostream>

#include "my_project/analytical.h"
#include "my_project/numerical.h"
#include "my_project/utils.h"

namespace acv {

void PrintAnalyticalAnalysis(const AnalyticalResult& res) {
    std::cout << "\n=== АНАЛИТИЧЕСКИЙ КРИТЕРИЙ УСТОЙЧИВОСТИ (Eq. 17) ===\n";

    if (res.is_stable) {
        std::cout << "[СТАТУС] СИСТЕМА УСТОЙЧИВА\n";
        std::cout << "Запас устойчивости (Margin): "
                  << std::abs(res.stability_margin) << "\n";
    } else {
        std::cout << "[СТАТУС] ВНИМАНИЕ: СИСТЕМА НЕУСТОЙЧИВА\n";
        std::cout << "Превышение порога: " << res.stability_margin << "\n";
    }

    std::cout << "--------------------------------------------------\n";
    std::cout << "Пневматический вклад (Term 1): " << res.pneumatic_term
              << "\n";
    std::cout << "Геометрический вклад (Term 2): " << res.geometric_term
              << "\n";

    if (res.influence_ratio > 1.0) {
        std::cout << "Доминирующий фактор: Характеристики нагнетателя и объем "
                     "подушки.\n";
    } else {
        std::cout
            << "Доминирующий фактор: Геометрия ограждения и угол наклона.\n";
    }

    if (!res.is_stable) {
        std::cout << "\nРЕКОМЕНДАЦИЯ: ";
        if (std::abs(res.pneumatic_term) < res.geometric_term) {
            std::cout << "Попробуйте увеличить рабочий угол phi или уменьшить "
                         "давление p0.\n";
        } else {
            std::cout << "Попробуйте увеличить объем подушки W0 или выбрать "
                         "нагнетатель с более крутой характеристикой.\n";
        }
    }
    std::cout << "==================================================\n";
}

void PrintResults(const StabilityResult& res) {
    std::cout << "--- Анализ устойчивости ---" << std::endl;
    std::cout << "Статус: " << (res.is_stable ? "УСТОЙЧИВА" : "НЕУСТОЙЧИВА")
              << std::endl;
    std::cout << "Макс. вещественная часть: " << res.max_real_part << std::endl;

    for (const auto& ev : res.eigenvalues) {
        std::cout << "Корень: " << ev.real() << " + " << ev.imag() << "i"
                  << std::endl;
    }

    std::cout << "\nКолебательные моды:" << std::endl;

    for (const auto& mode : res.oscillation_modes) {
        const auto& ev = mode.eigenvalue;

        std::cout << "  λ = " << ev.real() << " + " << ev.imag() << "i"
                  << std::endl;

        std::cout << "    Период T = " << mode.period << " c" << std::endl;

        std::cout << "    Лог. декремент χ = " << mode.logarithmic_decrement
                  << std::endl;

        std::cout << "    K = " << mode.decay_ratio << std::endl;

        std::cout << "    "
                  << (mode.logarithmic_decrement > 0 ? "Затухающие колебания"
                                                     : "Нарастающие колебания")
                  << std::endl;
    }

    std::cout << std::endl;
}

}  // namespace acv

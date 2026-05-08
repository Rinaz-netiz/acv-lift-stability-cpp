// src/utils.cpp
#include <iomanip>
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

void PrintVehicleTable(const VehicleData& data) {
    // Настройка формата вывода
    std::cout << std::fixed << std::setprecision(4);
    const int w = 30;  // Ширина колонки для имени переменной
    const std::string line(60, '-');

    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << " TABLE 1: INITIAL DATA AND PARAMETERS FOR " << data.name
              << "\n";
    std::cout << line << "\n";
    std::cout << std::left << std::setw(w) << "Property" << "Value" << "\n";
    std::cout << line << "\n";

    // 1. Геометрия и масса
    std::cout << std::setw(w) << "m (Mass half-ACV)" << data.geometry.m
              << " kg\n";
    std::cout << std::setw(w) << "L (Length)" << data.geometry.L << " m\n";
    std::cout << std::setw(w) << "l (Equivalent rod len)" << data.geometry.l
              << " m\n";
    std::cout << std::setw(w) << "S (Cushion area)" << data.geometry.S
              << " m2\n";
    std::cout << std::setw(w) << "W0 (Cushion volume)" << data.geometry.W0
              << " m3\n";

    // 2. Равновесные параметры
    std::cout << std::setw(w) << "p0 (Equilibrium pressure)" << data.eq.p0
              << " Pa\n";
    std::cout << std::setw(w) << "Q0 (Equilibrium flow)" << data.eq.Q0
              << " m3/s\n";
    std::cout << std::setw(w) << "phi0 (Seal angle)" << data.eq.phi0
              << " rad\n";

    // 3. Линеаризованные коэффициенты (как на фото)
    // dQ/dp в точке 0
    std::cout << std::setw(w) << "dQ/dp |_0" << data.der.dQin_dp
              << " m3/(s*Pa)\n";

    // Нормализованная жесткость ограждения (1/L * dM/dphi)
    double norm_dM_dphi =
        (data.geometry.L > 1e-6) ? (data.der.dM_dphi / data.geometry.L) : 0.0;
    std::cout << std::setw(w) << "1/L * dM/dphi |_0" << norm_dM_dphi << " N\n";

    // Демпфирование
    std::cout << std::setw(w) << "dM/dphi_dot |_0" << data.der.dM_dphidot
              << " N*m*s/rad\n";
    std::cout << std::setw(w) << "dY/dH_dot |_0" << data.der.dY_dHdot
              << " N*s/m\n";

    // Инерция
    std::cout << std::setw(w) << "I_phi (Seal inertia)" << data.geometry.I_phi
              << " kg*m2\n";

    std::cout << std::string(60, '=') << "\n" << std::endl;
}

}  // namespace acv

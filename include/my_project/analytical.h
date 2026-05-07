// include/my_project/analytical.h
#pragma once
#include "vehicle_data.h"

namespace acv {

struct AnalyticalResult {
    bool is_stable;
    double
        stability_margin;  // Запас устойчивости (насколько далеко мы от нуля)
    double pneumatic_term;  // Вклад нагнетателя и сжимаемости
    double geometric_term;  // Вклад геометрии ограждения
    double influence_ratio;  // Соотношение сил
};

AnalyticalResult VerifyAnalyticalDetailed(const VehicleData& data);

/**
 * Проверяет достаточное условие устойчивости (уравнение 32 из статьи).
 * Возвращает true, если система устойчива.
 */
bool AnalyticalVerification(const VehicleData& data);

}  // namespace acv

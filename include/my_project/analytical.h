// include/my_project/analytical.h
#pragma once
#include "vehicle_data.h"

namespace acv {

/**
 * Проверяет достаточное условие устойчивости (уравнение 32 из статьи).
 * Возвращает true, если система устойчива.
 */
bool analyticalVerification(const VehicleData& data);

}  // namespace acv

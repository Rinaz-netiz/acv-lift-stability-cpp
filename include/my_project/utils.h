// include/my_project/utils.h
#pragma once

#include "my_project/analytical.h"
#include "my_project/numerical.h"

namespace acv {

void PrintAnalyticalAnalysis(const AnalyticalResult& res);

void PrintResults(const StabilityResult& res);

void PrintVehicleTable(const VehicleData& data);

}  // namespace acv

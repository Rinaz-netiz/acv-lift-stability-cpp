#include <cmath>

#include "my_project/analytical.h"
#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {

AnalyticalResult VerifyAnalyticalDetailed(const VehicleData& data) {
    using namespace constants;

    // Слагаемое 1: Пневматическая составляющая (нагнетатель + давление)
    // Характеризует, как изменение давления "успокаивает" систему через расход
    double pneumatic = (data.l / 2.0) * (kN * kPAtm / data.W0) *
                       (data.dQin_dp - data.Q0 / (2.0 * data.p0));

    // Слагаемое 2: Геометрическая составляющая (эффект изменения зазора)
    // Характеризует, как поворот ограждения меняет площадь истечения
    double geometric =
        kChi * std::sqrt(2.0 * data.p0 / kRhoAir) * std::sin(data.phi0);

    AnalyticalResult res;
    res.stability_margin = pneumatic + geometric;  // Сумма должна быть < 0
    res.is_stable = res.stability_margin < 0;
    res.pneumatic_term = pneumatic;
    res.geometric_term = geometric;

    // Отношение: какой фактор доминирует
    res.influence_ratio = std::abs(pneumatic / geometric);

    return res;
}

bool AnalyticalVerification(const VehicleData& data) {
    using namespace constants;

    const double term1 = (data.l / 2.0) * (kN * kPAtm / data.W0) *
                         (data.dQin_dp - data.Q0 / (2.0 * data.p0));

    const double term2 =
        kChi * std::sqrt(2.0 * data.p0 / kRhoAir) * std::sin(data.phi0);

    return (term1 + term2) < 0;  // true → устойчиво
}

}  // namespace acv

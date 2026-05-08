// src/analytical.cpp
#include <cmath>

#include "my_project/analytical.h"
#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {

AnalyticalResult VerifyAnalyticalDetailed(const VehicleData& data) {
    using namespace constants;

    // Пневматическое слагаемое (используем data.der и data.eq)
    double pneumatic = (data.geometry.l / 2.0) *
                       (kN * kPAtm / data.geometry.W0) *
                       (data.der.dQin_dp - data.eq.Q0 / (2.0 * data.eq.p0));

    // Геометрическое слагаемое
    double geometric =
        kChi * std::sqrt(2.0 * data.eq.p0 / kRhoAir) * std::sin(data.eq.phi0);

    AnalyticalResult res;
    res.stability_margin = pneumatic + geometric;
    res.is_stable = res.stability_margin < 0;
    res.pneumatic_term = pneumatic;
    res.geometric_term = geometric;
    res.influence_ratio =
        (std::abs(geometric) > 1e-9) ? std::abs(pneumatic / geometric) : 0.0;

    return res;
}

bool AnalyticalVerification(const VehicleData& data) {
    using namespace constants;

    const double term1 = (data.geometry.l / 2.0) *
                         (kN * kPAtm / data.geometry.W0) *
                         (data.der.dQin_dp - data.eq.Q0 / (2.0 * data.eq.p0));

    const double term2 =
        kChi * std::sqrt(2.0 * data.eq.p0 / kRhoAir) * std::sin(data.eq.phi0);

    return (term1 + term2) < 0;
}

}  // namespace acv

#include "my_project/vehicle_data.h"
#include "my_project/analytical.h"
#include <cmath>

namespace acv {

bool analyticalVerification(const VehicleData& data) {
    const double term1 =
        (data.l / 2.0) *
        (data.n * data.p_atm / data.W0) *
        (data.dQin_dp - data.Q0 / (2.0 * data.p0));

    const double term2 =
        data.chi * std::sqrt(2.0 * data.p0 / data.rho_air) *
        std::sin(data.phi0);

    return (term1 + term2) < 0;  // true → устойчиво
}

} // namespace acv

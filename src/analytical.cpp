#include <cmath>

#include "my_project/analytical.h"
#include "my_project/constants.h"
#include "my_project/vehicle_data.h"

namespace acv {

bool analyticalVerification(const VehicleData& data) {
    using namespace constants;

    const double term1 = (data.l / 2.0) * (kN * kPAtm / data.W0) *
                         (data.dQin_dp - data.Q0 / (2.0 * data.p0));

    const double term2 =
        kChi * std::sqrt(2.0 * data.p0 / kRhoAir) * std::sin(data.phi0);

    return (term1 + term2) < 0;  // true → устойчиво
}

}  // namespace acv

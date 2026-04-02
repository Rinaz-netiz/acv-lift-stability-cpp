#include "my_project/vehicle_data.h"
#include "my_project/analytical.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <vehicle.json>\n";
        return 1;
    }

    try {
        auto v = acv::loadVehicleFromJson(argv[1]);
        bool stable = acv::analyticalVerification(v);
        std::cout << "Analytical stability check: "
                  << (stable ? "STABLE" : "UNSTABLE") << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

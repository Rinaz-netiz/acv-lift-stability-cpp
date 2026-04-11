// src/main.cpp
#include <iostream>

#include "my_project/analytical.h"
#include "my_project/numerical.h"

int main(int argc, char* argv[]) {
    auto v = acv::loadVehicleFromJson(argv[1]);

    // Аналитическая модель
    bool stable_analytical = acv::analyticalVerification(v);

    // Численная модель
    auto coeffs = acv::computeFullCharacteristicPolynomial(v);
    bool stable_numerical = acv::isStableNumerical(coeffs);

    std::cout << "Analytical: " << (stable_analytical ? "STABLE" : "UNSTABLE")
              << "\n";
    std::cout << "Numerical:  " << (stable_numerical ? "STABLE" : "UNSTABLE")
              << "\n";
}

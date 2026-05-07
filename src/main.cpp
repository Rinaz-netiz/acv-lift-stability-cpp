// src/main.cpp
#include <iostream>

#include "my_project/analytical.h"
#include "my_project/numerical.h"
#include "my_project/utils.h"

int main(int argc, char* argv[]) {
    auto v = acv::loadVehicleFromJson(argv[1]);

    try {
        auto analytical_res = acv::VerifyAnalyticalDetailed(v);
        acv::PrintAnalyticalAnalysis(analytical_res);
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
    }

    try {
        auto res = acv::AnalyzeStabilitySimple(v);
        std::cout << "Упрощенная модель:\n";
        acv::PrintResults(res);
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
    }

    try {
        auto res = acv::AnalyzeStabilityFull(v);
        std::cout << "Полная модель:\n";
        acv::PrintResults(res);
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
    }
}

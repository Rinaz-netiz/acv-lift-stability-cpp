cmake_minimum_required(VERSION 3.14)
project(ACV-LiftStability LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Используем vcpkg для nlohmann/json
find_package(nlohmann_json CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(acv-app
    src/main.cpp
    src/vehicle_data.cpp
    src/analytical.cpp
    src/numerical.cpp
    src/utils.cpp
)

target_include_directories(acv-app PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(acv-app PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(acv-app PRIVATE Eigen3::Eigen nlohmann_json::nlohmann_json)

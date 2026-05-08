cmake_minimum_required(VERSION 3.14)
project(ACV-LiftStability LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Опция для включения/выключения тестов
option(BUILD_TESTS "Build tests" ON)

# Используем vcpkg для зависимостей
find_package(nlohmann_json CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)

# ============================================
# Библиотека (чтобы переиспользовать код в тестах)
# ============================================
add_library(acv-lib STATIC
    src/vehicle_data.cpp
    src/analytical.cpp
    src/numerical.cpp
    src/utils.cpp
)

target_include_directories(acv-lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(acv-lib PUBLIC
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
)

# ============================================
# Основное приложение
# ============================================
add_executable(acv-app
    src/main.cpp
)

target_link_libraries(acv-app PRIVATE acv-lib)

# ============================================
# Тесты
# ============================================

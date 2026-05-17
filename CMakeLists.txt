cmake_minimum_required(VERSION 3.14)
project(ACV-LiftStability LANGUAGES CXX)

set(CMAKE_TOOLCHAIN_FILE "$ENV{HOME}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "vcpkg toolchain")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
if(BUILD_TESTS)
    enable_testing()

    # Найти GTest (через vcpkg или системный)
    find_package(GTest CONFIG REQUIRED)

    configure_file(
        ${CMAKE_SOURCE_DIR}/test_data/vehicle_a.json
        ${CMAKE_CURRENT_BINARY_DIR}/test_data/vehicle_a.json
        COPYONLY
    )

    # Исполняемый файл с тестами
    add_executable(acv-tests
        tests/test_vehicle_data.cpp
        tests/test_analytical.cpp
        tests/test_numerical.cpp
        tests/test_integration.cpp
    )

    target_link_libraries(acv-tests PRIVATE
        acv-lib
        GTest::gtest
        GTest::gtest_main
    )

    # Автоматическое обнаружение тестов
    include(GoogleTest)
    gtest_discover_tests(acv-tests)

    # Альтернативный способ (если gtest_discover_tests не работает)
    # add_test(NAME acv-tests COMMAND acv-tests)
endif()

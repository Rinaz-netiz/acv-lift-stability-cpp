# ACV Lift Stability Analyzer (C++)

Численный и аналитический анализ устойчивости подъемной системы судна на воздушной подушке (ACV/SES) с надувными бортовыми ограждениями.

Проект реализует:

- расчет равновесного состояния воздушной подушки;
- линеаризацию системы;
- анализ устойчивости по собственным значениям;
- расчет периодов и декрементов колебаний;
- аналитическую проверку достаточного условия устойчивости;
- загрузку и сохранение моделей КВП в JSON.

Модель основана на статье:

V. Shabarov, F. Peplin, P. Kalyasov, V. Shaposhnikov

"Analytical and numerical investigation of the lift system stability of the air cushion vehicle fitted with closed inflated side seals"

Applied Ocean Research, 2022

https://doi.org/10.1016/j.apor.2022.103045

---

# Возможности

## 1. Расчет равновесия

Автоматически вычисляются:

- давление в подушке `p0`
- расход воздуха `Q0`
- угол отклонения ограждения `phi0`
- зазор под ограждением `h_gap0`
- положение корпуса `H0`

На основе:

- характеристики нагнетателя
- расхода через щель
- момента ограждения

---

## 2. Линеаризация системы

Вычисляются производные:

- `dQin/dp`
- `dQout/dp`
- `dQout/dH`
- `dQout/dphi`
- `dM/dphi`
- `dM/dphidot`
- `dY/dHdot`

---

## 3. Упрощенная модель (3×3)

Состояния:

- Δp
- Δφ
- Δφ̇

Функция:

```cpp
acv::AnalyzeStabilitySimple(...)
```

---

## 4. Полная модель (5×5)

Состояния:

- Δp
- ΔH
- ΔḢ
- Δφ
- Δφ̇

Функция:

```cpp
acv::AnalyzeStabilityFull(...)
```

---

## 5. Аналитическая проверка устойчивости

Реализовано достаточное условие устойчивости из статьи.

Функции:

```cpp
acv::AnalyticalVerification(...)
acv::VerifyAnalyticalDetailed(...)
```

---

## 6. Анализ колебательных мод

Для каждой комплексной моды вычисляются:

- период колебаний
- логарифмический декремент
- коэффициент затухания

---

# Структура проекта

```text
acv-lift-stability-cpp/
│
├── data/
│   └── vehicle_A.json
│
├── include/my_project/
│   ├── analytical.h
│   ├── constants.h
│   ├── numerical.h
│   ├── utils.h
│   └── vehicle_data.h
│
├── src/
│   ├── analytical.cpp
│   ├── main.cpp
│   ├── numerical.cpp
│   ├── utils.cpp
│   └── vehicle_data.cpp
│
├── test_data/
│   └── vehicle_a.json
│
├── tests/
│   ├── test_analytical.cpp
│   ├── test_integration.cpp
│   ├── test_numerical.cpp
│   └── test_vehicle_data.cpp
│
├── CMakeLists.txt
└── README.md
```

---

# Требования

- C++17
- CMake ≥ 3.15
- Eigen3
- nlohmann/json
- GoogleTest (для тестов)

---

# Установка зависимостей

## Через vcpkg

```bash
vcpkg install eigen3
vcpkg install gtest
vcpkg install nlohmann-json
```

---

# Сборка проекта

## Сборка с тестами

```bash
mkdir build
cd build

cmake .. \
  -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake

cmake --build .
```

---

## Запуск тестов

```bash
ctest --output-on-failure
```

или напрямую:

```bash
./acv-tests
```

Запуск отдельного теста:

```bash
./acv-tests \
  --gtest_filter=VehicleDataTest.EquilibriumFlowBalance
```

---

# Запуск приложения

Из директории `build`:

```bash
./acv-app ../data/vehicle_A.json
```

---

# Формат входных данных

Входные параметры задаются в JSON.

Пример:

```json
{
  "name": "Vehicle A",

  "geometry": {
    "m": 2050.0,
    "L": 10.5,
    "l": 0.7,
    "S": 20.0,
    "W0": 14.0,
    "I_phi": 9.35,
    "H_init": 0.0
  },

  "blower": {
    "p_table": [0.0, 1000.0, 1005.525, 1010.0, 1500.0],
    "Q_table": [7.721, 5.72105, 5.71, 5.70105, 4.721]
  },

  "seal_moment": {
    "phi": [0.0, 0.1, 0.18, 0.19, 0.2, 0.3],
    "M_per_L": [0.0, -50.0, -224.5736, -246.3536, -268.1336, -485.9336]
  },

  "seal_damping_curve": {
    "h_gap_table": [0.0, 0.002, 0.006, 0.014, 0.1],
    "dM_dphidot_table": [-50.0, -30.0, -14.0, -5.0, -1.0]
  },

  "cushion_damping_curve": {
    "Sgap_over_S_table": [
      0.0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.02
    ],
    "D_table": [
      120000, 65000, 38000, 26000, 21000, 18500, 17000, 15000, 12000, 8000
    ]
  }
}
```

---

# Использование API

## Загрузка модели

```cpp
acv::VehicleData vehicle =
    acv::LoadVehicleFromJson("data/vehicle_A.json");
```

---

## Численный анализ устойчивости

### Упрощенная модель

```cpp
auto result =
    acv::AnalyzeStabilitySimple(vehicle);
```

### Полная модель

```cpp
auto result =
    acv::AnalyzeStabilityFull(vehicle);
```

---

## Аналитическая проверка

```cpp
auto analytical =
    acv::VerifyAnalyticalDetailed(vehicle);
```

---

# Интерпретация результатов

## Собственные значения

Система устойчива если:

```text
Re(λ) < 0
```

для всех собственных значений.

---

## Колебательная мода

Для комплексного корня:

```text
λ = σ + iω
```

вычисляются:

Период:

```text
T = 2π / |ω|
```

Логарифмический декремент:

```text
χ = -σT
```

Коэффициент затухания:

```text
K = exp(χ)
```

---

## Интерпретация декремента

| Условие | Поведение          |
| ------- | ------------------ |
| χ > 0   | колебания затухают |
| χ < 0   | колебания растут   |
| K > 1   | устойчивая мода    |
| K < 1   | неустойчивая мода  |

---

# Основные классы

## VehicleData

Полное описание КВП:

```cpp
acv::VehicleData
```

Содержит:

- геометрию
- характеристики нагнетателя
- характеристики ограждения
- равновесное состояние
- производные линеаризации

---

## StabilityResult

Результат численного анализа:

```cpp
acv::StabilityResult
```

Содержит:

- собственные значения
- максимальную действительную часть
- список колебательных мод
- флаг устойчивости

---

## AnalyticalResult

Результат аналитической проверки:

```cpp
acv::AnalyticalResult
```

Содержит:

- флаг устойчивости
- запас устойчивости
- пневматическое слагаемое
- геометрическое слагаемое

---

# Используемые методы

В проекте используются:

- линейная интерполяция
- бисекция для поиска равновесия
- линеаризация системы
- анализ собственных значений
- критерий устойчивости
- динамический анализ колебаний

---

# Тестирование

Покрытие включает:

- интерполяцию
- вычисление равновесия
- проверку производных
- аналитическую устойчивость
- численный анализ
- интеграционные тесты
- граничные случаи

---

# Лицензия

MIT License

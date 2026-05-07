# ACV Lift Stability Analyzer (C++)

Численный и аналитический анализ устойчивости подъемной системы судна на воздушной подушке (ACV/SES) с надувными бортовыми ограждениями.

Проект реализует:

- линейный анализ устойчивости;
- расчет собственных значений системы;
- вычисление периодов и декрементов колебаний;
- аналитическую проверку условия устойчивости из статьи:
  "Analytical and numerical investigation of the lift system stability of the air cushion vehicle fitted with closed inflated side seals".

---

## Возможности

---

Проект поддерживает:

1. Упрощенную модель (3×3)

Состояния:

- Δp
- Δφ
- Δφ̇

Функция:

```cpp
AnalyzeStabilitySimple(...)
```

2. Полную модель (5×5)

Состояния:

- Δp
- ΔH
- ΔḢ
- Δφ
- Δφ̇

Функция:

```cpp
AnalyzeStabilityFull(...)
```

3. Аналитическую проверку устойчивости

Условие (32) из статьи.

Функции:

```cpp
AnalyticalVerification(...)
VerifyAnalyticalDetailed(...)
```

4. Расчет колебательных характеристик

Для каждой комплексной моды вычисляются:

- период колебаний;
- логарифмический декремент;
- коэффициент затухания.

---

## Структура проекта

---

```text
acv-lift-stability-cpp/
│
├── data/
│   ├── vehicle_A.json
│   ├── vehicle_B.json
│   └── vehicle_C.json
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
├── tests/
│   ├── test_analytical.cpp
│   ├── test_integration.cpp
│   ├── test_numerical.cpp
│   └── test_vehicle_data.cpp
│
├── CMakeLists.txt
├── compile_commands.json
└── README.md
```

---

## Требования

---

- C++17
- CMake ≥ 3.15
- Eigen3
- GoogleTest (для тестов)
- vcpkg (рекомендуется)

---

## Установка зависимостей

**Через vcpkg**

Установить:

```bash
vcpkg install eigen3
vcpkg install gtest
```

---

## Сборка проекта

---

```bash
# Сборка с тестами (по умолчанию)
mkdir build && cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build .

# Запуск всех тестов
ctest --output-on-failure

# Или напрямую
./acv-tests

# Запуск конкретного теста
./acv-tests --gtest_filter=VehicleDataTest.InitComputesDerivativesCorrectly

# Сборка БЕЗ тестов
cmake .. -DBUILD_TESTS=OFF
cmake --build .
```

После сборки появится исполняемый файл.

---

## Запуск

---

Пример:

```bash
./acv-app ../data/vehicle_A.json
```

---

## Формат входных данных

---

Параметры судна задаются в JSON.

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
    "I_phi": 9.35
  },

  "equilibrium": {
    "p0": 1005.5,
    "Q0": 5.71,
    "phi0": 0.19
  },

  "derivatives": {
    "dQin_dp": -0.002,
    "dM_dphi_per_L": -2178.0,
    "dM_dphidot": -3.0,
    "dY_dHdot_factor": -17836.0
  }
}
```

---

## Параметры VehicleData

---

### Геометрия

| Поле  | Описание                       | Единицы |
| ----- | ------------------------------ | ------- |
| m     | масса                          | kg      |
| L     | длина                          | m       |
| l     | эквивалентная длина ограждения | m       |
| S     | площадь воздушной подушки      | m²      |
| W0    | объем воздушной подушки        | m³      |
| I_phi | момент инерции ограждения      | kg·m²   |

### Равновесное состояние

| Поле | Описание                   |
| ---- | -------------------------- |
| p0   | давление в подушке         |
| Q0   | расход воздуха             |
| phi0 | угол отклонения ограждения |

### Производные

| Поле            | Описание                               |
| --------------- | -------------------------------------- |
| dQin_dp         | производная характеристики нагнетателя |
| dM_dphi_per_L   | жесткость ограждения                   |
| dM_dphidot      | демпфирование ограждения               |
| dY_dHdot_factor | коэффициент демпфирования подушки      |

---

## Использование API

---

### Загрузка данных

```cpp
auto vehicle = acv::loadVehicleFromJson("data/vehicle_A.json");
```

---

### Численный анализ

```cpp
auto result = acv::AnalyzeStabilityFull(vehicle);
```

или

```cpp
auto result = acv::AnalyzeStabilitySimple(vehicle);
```

---

### Аналитическая проверка

```cpp
auto analytical =
    acv::VerifyAnalyticalDetailed(vehicle);
```

---

## Интерпретация результатов

---

### Собственные значения

Система устойчива если:

```text
Re(λ) < 0
```

для всех корней.

---

### Колебательные моды

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

### Интерпретация декремента

| Условие | Поведение          |
| ------- | ------------------ |
| χ > 0   | колебания затухают |
| χ < 0   | колебания растут   |
| K > 1   | устойчивая мода    |
| K < 1   | неустойчивая мода  |

---

---

## Используемая теория

---

В основе проекта:

- линеаризация системы;
- матрица Якоби;
- собственные значения линейной системы;
- критерий Рауса–Гурвица;
- анализ затухающих колебаний.

---

## Источник модели

---

V. Shabarov, F. Peplin, P. Kalyasov, V. Shaposhnikov

"Analytical and numerical investigation of the lift system stability of the air cushion vehicle fitted with closed inflated side seals"

Applied Ocean Research, 2022.

DOI:
https://doi.org/10.1016/j.apor.2022.103045

---

## Лицензия

---

MIT License

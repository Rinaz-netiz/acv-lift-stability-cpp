# ACV Lift Stability Analysis

Численное и аналитическое исследование устойчивости системы воздушной подушки
судна на воздушной подушке с надувными бортовыми уплотнениями.

## Математическая модель

Основано на статье:
Analytical and numerical investigation of the lift system stability
of the air cushion vehicle fitted with closed inflated side seals.
Vasiliy Shabarov, Fedor Peplin, Pavel Kalyasov, Vitaliy Shaposhnikov (2022), Applied Ocean Research

### Характеристическое уравнение

`λ⁵ + a₄λ⁴ + a₃λ³ + a₂λ² + a₁λ + a₀ = 0`

### Критерий устойчивости (упрощённая модель)

`l/2 · (np_atm/W₀)(dQin/dp - Q₀/2p₀) + χ√(2p₀/ρ)sin(φ₀) < 0`

## Использование

```bash
mkdir build && cd build
cmake ..
make
./acv_stability ../data/vehicle_A.json
```

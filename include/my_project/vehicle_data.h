// include/my_project/vehicle_data.h
#pragma once
#include <string>
#include <vector>

namespace acv {

// Входные данные, задаются пользователем
struct VehicleGeometry {
    double m;  // Масса половины КВП, kg
    double L;  // Длина, m
    double l;  // Длина эквивалентного стержня (ограждения), m
    double S;  // Площадь воздушной подушки, m²
    double W0;      // Объём воздушной подушки, m³
    double I_phi;   // Момент инерции ограждения, kg·m²
    double H_init;  // Высота ЦТ в начальном состоянии (без подушки), m
};

// Характеристика нагнетателя Q_in(p) — Eq. (12)
// Задаётся как таблица точек (p_i, Q_i) — интерполируется линейно
// ИЛИ как линейная модель: Q = Q_max - k_fan * p
struct BlowerCurve {
    // Линейная модель
    double Q_max = 0.0;  // m³/s при p=0
    double k_fan = 0.0;  // m³/(s·Pa),  dQ/dp = -k_fan < 0

    // Табличная модель (если не пустая — используется вместо линейной)
    std::vector<double> p_table;  // Pa
    std::vector<double> Q_table;  // m³/s

    // Вычислить Q(p) и dQ/dp
    double Flow(double p) const;
    double Derivative(double p) const;  // dQ/dp < 0
};

// Момент сопротивления ограждения M(φ) — Fig. 13, Section 2.2.1
// Задаётся как таблица пар (φ_i, M_i) — интерполируется линейно
// M < 0 при φ > 0 (возвращающий момент)
struct SealMomentCurve {
    std::vector<double> phi_table;  // rad
    std::vector<double> M_per_L_table;  // N (момент на единицу длины L)

    // Вычислить M(φ)*L и dM/dφ в точке φ
    double Moment(double phi, double L) const;      // N·m
    double Derivative(double phi, double L) const;  // N·m/rad,  < 0
};

struct SealDampingCurve {
    std::vector<double> h_gap_table;       // m
    std::vector<double> dM_dphidot_table;  // N·m·s/rad (<0)

    double Value(double h_gap) const;  // линейная интерполяция
};

struct CushionDampingCurve {
    std::vector<double> Sgap_over_S_table;  // (-)
    std::vector<double> D_table;            // (-), как в Eq.(13)

    double Value(double Sgap_over_S) const;  // линейная интерполяция
};

// Результаты вычислений равновесия
struct EquilibriumState {
    double p0;    // Избыточное давление в подушке, Pa
    double phi0;  // Угол отклонения ограждения, rad
    double H0;    // Высота ЦТ при равновесии, m
    double Q0;    // Расход в равновесии, m³/s
    double h_gap0;  // Зазор под ограждением, m
};

struct LinearizedDerivatives {
    // Eq. (25)
    double dQout_dp;    // = Q0 / (2*p0),  > 0
    double dQout_dH;    // = χ*sqrt(2p0/ρ)*L,  > 0
    double dQout_dphi;  // = χ*sqrt(2p0/ρ)*L*l*sin(φ0),  > 0 при φ0>0

    double dQin_dp;  // < 0 (из характеристики нагнетателя)

    // Из Section 2.2.1 — производная момента ограждения по углу
    double dM_dphi;  // < 0

    // Из Section 2.2.2 — демпфирование ограждения (CFD)
    double dM_dphidot;  // < 0

    // Из Eq. (13) — демпфирование воздушной подушки (CFD)
    // dY/dH_dot = D(S_gap) * rho * (-Q_in/S) * S = -D*rho*Q_in
    double dY_dHdot;  // < 0
};

struct VehicleData {
    std::string name;

    // Входные данные
    VehicleGeometry geometry{};
    BlowerCurve blower{};
    SealMomentCurve seal_moment{};
    SealDampingCurve seal_damping_curve{};
    CushionDampingCurve cushion_damping_curve{};

    // Демпфирование ограждения dM/dφ̇ — из CFD (Section 2.2.2)
    // Задаётся как константа (в статье — рис. 9, значения из Table 1)
    double dM_dphidot_input = -1.0;  // N·m·s/rad,  < 0

    // Демпфирование подушки dY/dḢ — из Eq. (13) и CFD (рис. 4)
    // dY/dHdot = D(h_gap) * rho * (-Q_in/S) * S
    // D берётся из рис. 4 для h_gap0
    // Задаём как: dY_dHdot = D_coeff * rho * Q0  (с нужным знаком)
    double D_damping = 0.0;  // коэффициент D из рис. 4 (положительный)

    // Вычисляемые результаты
    EquilibriumState eq{};        // Равновесие
    LinearizedDerivatives der{};  // Производные

    void Init();  // Вычислить eq и der, проверить

   private:
    void ComputeEquilibrium();
    void ComputeDerivatives();
    void Validate() const;
    void SortingInputedCurve(std::vector<double>& arguments,
                             std::vector<double>& values);
};

// Загрузка / сохранение
VehicleData LoadVehicleFromJson(const std::string& filename);
void SaveVehicleToJson(const VehicleData& v, const std::string& filename);

}  // namespace acv

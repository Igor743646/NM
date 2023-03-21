#pragma once
#include <cmath>
#include <functional>
#include <unordered_map>
#include "utils.h"

class FDESbase {

protected:

    double L, R;            // границы отрезка сетки по x координате
    double T;               // граница отрезка времени по t координате
    double alpha, gamma;    // степени производных по x и t координатах соответственно
    double h, tau;          // шаги по сетки по x и t координатах соответственно
    double beta;            // коэффициент доли лево/правосторонней производных [-1; +1]

    std::function<double(double, double)> D;    // коэффициент диффузии при дробной производной по пространству
    std::function<double(double, double)> V;    // коэффициент сноса при производной первой степени
    std::function<double(double)> psi;          // начальное условие при t = 0, u(x, 0) = psi(x)
    std::function<double(double, double)> f;    // функция источник
    std::function<double(double)> phiL;         // граничное условие u(L, t) = phiL(t)
    std::function<double(double)> phiR;         // граничное условие u(R, t) = phiR(t)

    ull n, k;       // количество ячеек по x и t координатах соответственно

    std::unordered_map<ull, double> memo_g_alpha;
    std::unordered_map<ull, double> memo_g_gamma;

public:

    std::vector<std::vector<double>> result; // сетка с результатом моделирования/вычислений

    FDESbase(double, double, double, double, double, double, double, double,
        std::function<double(double, double)>, std::function<double(double, double)>, std::function<double(double)>,
        std::function<double(double, double)>, std::function<double(double)>, std::function<double(double)>);

    FDESbase(double, double, double, double, double, double, double, double,
        std::function<double(double, double)>, std::function<double(double, double)>, std::function<double(double)>,
        std::function<double(double, double)>);

    FDESbase(const FDESbase& p);

    void Info();

    double g(double, ull);

    double x_i(ull);
    double t_k(ull);

    double a(double, double);
    double b(double, double);
    double c(double, double);

    friend std::ostream& operator<<(std::ostream&, const FDESbase&);

};
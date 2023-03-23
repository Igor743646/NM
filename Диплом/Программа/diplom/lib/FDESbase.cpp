#include "FDESbase.h"

FDESbase::FDESbase(double _L, double _R, double _T, double _alpha, double _gamma, double _h, double _tau, double _beta,
    std::function<double(double, double)> _D, 
    std::function<double(double, double)> _V, 
    std::function<double(double)> _psi,
    std::function<double(double, double)> _f, 
    std::function<double(double)> _phiL,
    std::function<double(double)> _phiR,
    double _alpha_l, double _beta_l,
    double _alpha_r, double _beta_r) 
    : L(_L), R(_R), T(_T), alpha(_alpha), gamma(_gamma), h(_h), tau(_tau), 
    beta(_beta), D(_D), V(_V), psi(_psi), f(_f), phiL(_phiL), phiR(_phiR),
    alpha_l(_alpha_l), beta_l(_beta_l), alpha_r(_alpha_r), beta_r(_beta_r)
{
    n = (ull)((R - L) / h);
    k = (ull)(T / tau);
    result = std::vector<std::vector<double>>(k + 1, std::vector<double>(n + 1, 0.0));
    memo_g_alpha[0] = 1.0;
    memo_g_gamma[0] = 1.0;
}

FDESbase::FDESbase(const FDESbase& p) 
: L(p.L), R(p.R), T(p.T), alpha(p.alpha), gamma(p.gamma), h(p.h), tau(p.tau), 
beta(p.beta), D(p.D), V(p.V), psi(p.psi), f(p.f), phiL(p.phiL), phiR(p.phiR), 
alpha_l(p.alpha_l), beta_l(p.beta_l), alpha_r(p.alpha_r), beta_r(p.beta_r), n(p.n), k(p.k)
{
    result = p.result;
    memo_g_alpha[0] = 1.0;
    memo_g_gamma[0] = 1.0;
}

void FDESbase::Info() {

    Output::print("Info:");
    Output::print("\tx : [", L, "; ", R, "], t : [", 0, "; ", T, "]");
    Output::print("\th = ", h, "; tau = ", tau);
    Output::print("\tN = ", n);
    Output::print("\tT = ", k);
    Output::print("\tD * (tau^gamma) / (h^alpha) = ", D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha));
    Output::print("\tgamma / alpha = ", gamma / alpha);

    if (D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha) > gamma / alpha) {
        Output::print("May be problem with condition");
        Output::print(D(0, 0) * std::pow(tau, gamma) / std::pow(h, alpha), " ", gamma / alpha);
        Output::print("You may make tau <= ", std::pow(gamma * std::pow(h, alpha) / alpha / D(0, 0), 1.0/gamma));
        Output::print("Or make h >= ", std::pow(std::pow(tau, gamma) * alpha * D(0, 0) / gamma, 1.0/alpha));
        throw "ERROR: Bad conditions\n";
    }
    
    Output::print("\n");
}

double FDESbase::g(double a, ull j) {

    if (a != alpha && a != gamma) {
        throw "ERROR: g function using memoization only for \"alpha\" and \"gamma\". If you need more use comments below\n";

        /*
            if (std::abs(a) - (ull)std::abs(a) < 0.000001 && j >= a + 1) {
                return 0.0;
            }
            return (j % 2 == 0 ? 1 : -1) / (a + 1) / BETA((double)j + 1.0, a - (double)j + 1.0);
        */
    }

    if (a == alpha) {
        if (memo_g_alpha.find(j) != memo_g_alpha.end()) {
            return memo_g_alpha[j];
        } else {
            memo_g_alpha[j] = (j - 1.0 - alpha) / j * FDESbase::g(a, j - 1);
            return memo_g_alpha[j];
        }
    } 
    
    if (memo_g_gamma.find(j) == memo_g_gamma.end()) {
        memo_g_gamma[j] = (j - 1.0 - gamma) / j * FDESbase::g(a, j - 1);
    }

    return memo_g_gamma[j];
}

double FDESbase::x_i(ull i) {
    return L + i * h;
}

double FDESbase::t_k(ull k) {
    return k * tau;
}

double FDESbase::a(double x, double t) {
    return (1.0 + beta) * (D(x, t) / 2.0) * (std::pow(tau, gamma) / std::pow(h, alpha));
}

double FDESbase::b(double x, double t) {
    return (1.0 - beta) * (D(x, t) / 2.0) * (std::pow(tau, gamma) / std::pow(h, alpha));
}

double FDESbase::c(double x, double t) {
    return V(x, t) * std::pow(tau, gamma) / 2 / h;
}

std::ostream& operator<<(std::ostream& out, const FDESbase& p) {
    out << p.k + 1 << std::endl << p.h << std::endl << p.tau << std::endl << p.L << std::endl;
    for (ull t = 0; t < p.result.size(); t++) {
        for (ull i = 0; i < p.result[t].size(); i++) {
            out << p.result[t][i] << " ";
        }
        if (t < p.result.size() - 1)
            out << std::endl;
    }
    return out;
}
#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <limits>

using ull = unsigned long long;

namespace std {
    constexpr double PI = 3.14159265359;
}

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
    for (auto& a : v) {
        out << a << " ";
    }
    return out;
}

namespace Output {

    template<class T>
    void print(T object, bool end = true) {
        if (end == true)
            std::cout << object << std::endl;
        else
            std::cout << object;
    }

    template<class T, class... Args>
    void print(T&& first, Args... next) {
        print(first, false);
        print(next...);
    }

}

inline double GAMMA(double x) {
    return x > 0 ? std::tgamma(x) : -M_PI / (x * std::sin(M_PI * x) * std::tgamma(-x));
}

inline double BETA(double x, double y) {
    double res = GAMMA(x) * GAMMA(y) / GAMMA(x + y);
    if (std::isnan(res)) {
        // std::cout << x << " " << y << std::endl;
        // throw "ERROR: NAN value of BETA function\n";
        res = std::numeric_limits<double>::infinity();
    }
    return res;
}

template<class T>
std::vector<T> prefix_sum(const std::vector<T>& vec) {
    std::vector<T> result(vec.size(), vec[0]);

    for (ull i = 1; i < vec.size(); i++) {
        result[i] = vec[i] + result[i - 1];
    }

    return result;
}

template<class T>
ull bin_search(const std::vector<T>& v, T k) {
        ull l = 0, r = v.size() - 1, mid;

    while (r - l > 1) {

        mid = (r + l) / 2;

        if (k < v[mid]) {
            r = mid;
        }
        else if (k > v[mid]) {
            l = mid;
        }
        else {
            l = mid;
            r = mid;
        }

    }

    return k <= v[l] ? l : r;
}
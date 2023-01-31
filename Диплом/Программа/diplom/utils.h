#pragma once
#include <vector>
#include <iostream>
#include <cmath>



void print() {
    std::cout << std::endl;
}

template<class T>
void print(T obj) {
    std::cout << obj << std::endl;
}

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
    for (auto& a : v) {
        out << a << " ";
    }
    return out;
}

inline double GAMMA(double x) {
    return x > 0 ? std::tgamma(x) : -M_PI / (x * std::sin(M_PI * x) * std::tgamma(-x));
}

inline double BETA(double x, double y) {
    return GAMMA(x) * GAMMA(y) / GAMMA(x + y);
}
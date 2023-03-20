#pragma once
#include <vector>
#include <iostream>

namespace {
    constexpr double MECH_EPS = 0.0000001;
    using ull = unsigned long long;
}

struct Matrix {

    ull n;
    double** matrix;

    Matrix(ull _n) {
        n = _n;
        matrix = new double* [n];
        for (ull i = 0; i < n; i++) {
            matrix[i] = new double[n]{0};
        }
    }

    ~Matrix() {
        for (ull i = 0; i < n; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;
    }

    Matrix(ull _n, std::initializer_list<double> list) : Matrix(_n) {
        if (list.size() != _n * _n) {
            throw "Error: size of matrix not equal list size\n";
        }

        auto it = list.begin();

        for (ull i = 0; i < size(); i++) {
            for (ull j = 0; j < size(); j++) {
                matrix[i][j] = *it;
                it++;
            }
        }
    }

    Matrix(const Matrix& m) : Matrix(m.n) {
        for (ull i = 0; i < size(); i++) {
            for (ull j = 0; j < size(); j++) {
                matrix[i][j] = m[i][j];
            }
        }
    }

    ull size() const {
        return n;
    }

    void SwapLines(ull i, ull j) {
        if (i >= n || j >= n) {
            throw "Error: bad index in SwapLines\n";
        }

        std::swap(matrix[i], matrix[j]);
    }

    void SwapColumns(ull i, ull j) {
        if (i >= n || j >= n) {
            throw "Error: bad index in SwapColumns\n";
        }

        for (ull k = 0; k < size(); k++) {
            std::swap(matrix[k][i], matrix[k][j]);
        }
    }

    static Matrix E(ull _n) {
        Matrix result(_n);
        for (ull i = 0; i < _n; i++) {
            result[i][i] = 1;
        }
        return result;
    }

    std::vector<Matrix> LUFactorizing(int* count = nullptr) {
        Matrix P = E(size());
        Matrix L = E(size());
        Matrix U(*this);

        /*std::cout << P << std::endl;
        std::cout << L << std::endl;
        std::cout << U << std::endl;*/

        for (ull i = 0; i < size(); i++) {

            // 1. Находим строку с максимальным по модулю элементом.
            {
                ull k = i;
                double max = std::abs(U[i][i]);

                for (ull j = i + 1; j < size(); j++) {
                    if (std::abs(U[j][i]) > max) {
                        max = std::abs(U[j][i]);
                        k = j;
                    }
                }

                if (U[k][i] == 0) {
                    continue;
                }

                // 2. Меняем строки в U и обновляем L.
                if (k != i) {
                    P.SwapColumns(i, k);
                    L.SwapLines(i, k);
                    L.SwapColumns(i, k);
                    U.SwapLines(i, k);
                    if (count != nullptr) (*count) += 1;
                }
            }

            // 3. Алгоритм Гаусса
            for (ull j = i + 1; j < size(); j++) {
                double koef = U[j][i] / U[i][i];

                U[j][i] = 0;
                L[j][i] = koef;

                for (ull t = i + 1; t < size(); t++) {
                    U[j][t] -= koef * U[i][t];
                }
            }
        }

        // std::cout << P << std::endl;
        // std::cout << L << std::endl;
        // std::cout << U << std::endl << std::endl;

        return std::vector<Matrix>({ P, L, U });
    }

    std::vector<double> Solve(const std::vector<double>& b) {
        // A * x = b => P * L * U * x = b => L * U * x = P^(-1) * b = P^(T) * b

        if (b.size() != size()) {
            throw "Error: matrix size and vector size are not equal\n";
        }

        // 1. Делаем LU - разложение
        auto plu = LUFactorizing();

        // 2. Вычисляем P^(T) * b = b * P = y
        auto y = b * plu[0];

        // 3. Вычисляем L * z = y;
        std::vector<double> z(size(), 0.0);

        for (ull i = 0; i < size(); i++) {
            z[i] = y[i];

            for (ull j = 0; j < i; j++) {
                z[i] -= plu[1][i][j] * z[j];
            }

            z[i] /= plu[1][i][i];
        }

        // 4. Вычисляем U * x = z
        std::vector<double> x(size(), 0.0);

        for (long i = size() - 1; i >= 0; i--) {
            x[i] = z[i];

            for (long j = i + 1; j < size(); j++) {
                x[i] -= plu[2][i][j] * x[j];
            }

            x[i] /= plu[2][i][i];
        }

        return x;
    }

    double Determinant() {
        int count = 0;
        auto p = LUFactorizing(&count);
        double result = 1;

        for (ull i = 0; i < size(); i++) {
            result *= p[2][i][i];
        }

        return count % 2 == 0 ? result : -result;
    }

    Matrix Reverse() {
        Matrix result(size());

        std::vector<double> b(size(), 0.0);

        for (ull i = 0; i < size(); i++) {
            b[i] = 1;
            auto res = Solve(b);
            for (ull j = 0; j < size(); j++) {
                result[j][i] = res[j];
            }
            b[i] = 0;
        }

        return result;
    }

    // Matrix<T>& operator= (std::initializer_list<T> list) {
    //     if (list.size() != size() * size()) {
    //         throw "error";
    //     }

    //     auto it = list.begin();

    //     for (ull i = 0; i < size(); i++) {
    //         for (ull j = 0; j < size(); j++) {
    //             matrix[i][j] = *it;
    //             it++;
    //         }
    //     }

    //     return *this;
    // }

    Matrix& operator= (Matrix m) {
        if (m.size() * m.size() != size() * size()) {
            throw "Error: sizes are not equal\n";
        }

        for (ull i = 0; i < size(); i++) {
            for (ull j = 0; j < size(); j++) {
                matrix[i][j] = m[i][j];
            }
        }

        return *this;
    }

    friend Matrix operator*(const Matrix& m1, const Matrix& m2) {
        Matrix result(m1.size());

        for (ull i = 0; i < m1.size(); i++) {
            for (ull j = 0; j < m1.size(); j++) {
                for (ull k = 0; k < m1.size(); k++) {
                    result[i][j] += 
                        m1[i][k] * 
                        m2[k][j];
                }
            }
        }

        return result;
    }

    // friend std::vector<T> operator*(const Matrix<T>& m1, const std::vector<T>& m2) {
    //     if (m1.size() != m2.size()) {
    //         throw "bad thing";
    //     }

    //     std::vector<T> result(m1.size(), T());
    //     for (ull i = 0; i < m1.size(); i++) {
    //         for (ull j = 0; j < m1.size(); j++) {
    //             result[i] += m1[i][j] * m2[j];
    //         }
    //     }

    //     return result;
    // }

    friend std::vector<double> operator*(const std::vector<double>& m2, const Matrix& m1) {
        if (m1.size() != m2.size()) {
            throw "Error: matrix size and vector size are not equal\n";
        }

        std::vector<double> result(m1.size(), 0.0);
        for (ull i = 0; i < m1.size(); i++) {
            for (ull j = 0; j < m1.size(); j++) {
                result[i] += m1[j][i] * m2[j];
            }
        }

        return result;
    }

    double* operator[](const ull i) const {
        return matrix[i];
    }

    friend std::ostream& operator<< (std::ostream& out, const Matrix& m) {
        for (ull i = 0; i < m.size(); i++) {
            for (ull j = 0; j < m.size(); j++) {
                out << m[i][j] << ' ';
            }
            out << std::endl;
        }
        return out;
    }

};


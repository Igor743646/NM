#include <iostream>
#include <initializer_list>
#include <vector>

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

template <class T>
class Matrix {

    static constexpr double MECH_EPS = 0.0000001;
    std::vector<std::vector<T>> matrix;
    std::vector<Matrix<T>> plu;

public:

    Matrix(size_t _n) {
        matrix = std::vector<std::vector<T>>(_n, std::vector<T>(_n, T()));
    }

    Matrix(size_t _n, std::vector<std::vector<T>>& matr) : Matrix(_n) {
        for (size_t i = 0; i < _n; i++) {
            for (size_t j = 0; j < _n; j++) {
                matrix[i][j] = matr[i][j];
            }
        }
    }

    Matrix(size_t _n, std::initializer_list<T> list) : Matrix(_n) {
        if (list.size() != _n * _n) {
            throw "error";
        }

        auto it = list.begin();

        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                matrix[i][j] = *it;
                it++;
            }
        }
    }

    size_t size() const {
        return matrix.size();
    }

    void SwapLines(size_t i, size_t j) {
        for (size_t k = 0; k < size(); k++) {
            std::swap(matrix[i][k], matrix[j][k]);
        }
    }

    void SwapColumns(size_t i, size_t j) {
        for (size_t k = 0; k < size(); k++) {
            std::swap(matrix[k][i], matrix[k][j]);
        }
    }

    Matrix<T> E(size_t _n) {
        Matrix<T> result(_n);
        for (size_t i = 0; i < _n; i++) {
            result[i][i] = 1;
        }
        return result;
    }

    std::vector<Matrix<T>> LUFactorizing(int* count = nullptr) {
        Matrix<T> P = E(size());
        Matrix<T> L = E(size());
        Matrix<T> U(size(), matrix);

        /*std::cout << P << std::endl;
        std::cout << L << std::endl;
        std::cout << U << std::endl;*/

        for (size_t i = 0; i < size(); i++) {

            // 1. Находим строку с максимальным по модулю элементом.
            {
                size_t k = i;
                T max = std::abs(U[i][i]);

                for (size_t j = i + 1; j < size(); j++) {
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
            for (size_t j = i + 1; j < size(); j++) {
                double koef = U[j][i] / U[i][i];

                U[j][i] = 0;
                L[j][i] = koef;

                for (size_t t = i + 1; t < size(); t++) {
                    U[j][t] -= koef * U[i][t];
                }
            }
        }

        /*std::cout << P << std::endl;
        std::cout << L << std::endl;
        std::cout << U << std::endl;*/

        return std::vector<Matrix<T>>({ P, L, U });
    }

    std::vector<T> Solve(const std::vector<T>& b) {
        // A * x = b => P * L * U * x = b => L * U * x = P^(-1) * b = P^(T) * b

        if (b.size() != size()) throw "размерность не совпадает";

        // 1. Делаем LU - разложение
        if (plu.size() == 0) {
            plu = LUFactorizing();
        }

        // 2. Вычисляем P^(T) * b = b * P = y
        auto y = b * plu[0];

        // 3. Вычисляем L * z = y;
        std::vector<T> z(size(), T());

        for (size_t i = 0; i < size(); i++) {
            z[i] = y[i];

            for (size_t j = 0; j < i; j++) {
                z[i] -= plu[1][i][j] * z[j];
            }

            z[i] /= plu[1][i][i];
        }

        // 4. Вычисляем U * x = z
        std::vector<T> x(size(), T());

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

        for (size_t i = 0; i < size(); i++) {
            result *= p[2][i][i];
        }

        return count % 2 == 0 ? result : -result;
    }

    Matrix<T> Reverse() {
        Matrix<T> result(size());

        std::vector<T> b(size(), 0);

        for (size_t i = 0; i < size(); i++) {
            b[i] = 1;
            auto res = Solve(b);
            for (size_t j = 0; j < size(); j++) {
                result[j][i] = res[j];
            }
            b[i] = 0;
        }

        return result;
    }

    Matrix<T>& operator= (std::initializer_list<T> list) {
        if (list.size() != size() * size()) {
            throw "error";
        }

        auto it = list.begin();

        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                matrix[i][j] = *it;
                it++;
            }
        }

        return *this;
    }
    friend Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2) {
        std::vector<std::vector<T>> vec(m1.size(), std::vector<T>(m1.size(), T()));
        for (size_t i = 0; i < m1.size(); i++) {
            for (size_t j = 0; j < m1.size(); j++) {
                for (size_t k = 0; k < m1.size(); k++) {
                    vec[i][j] += 
                        m1[i][k] * 
                        m2[k][j];
                }

                if (std::abs(vec[i][j]) < 2 * MECH_EPS) vec[i][j] = 0;
            }
        }

        return Matrix<T>(m1.size(), vec);
    }
    friend std::vector<T> operator*(const Matrix<T>& m1, const std::vector<T>& m2) {
        if (m1.size() != m2.size()) {
            throw "bad thing";
        }

        std::vector<T> result(m1.size(), T());
        for (size_t i = 0; i < m1.size(); i++) {
            for (size_t j = 0; j < m1.size(); j++) {
                result[i] += m1[i][j] * m2[j];
            }
        }

        return result;
    }
    friend std::vector<T> operator*(const std::vector<T>& m2, const Matrix<T>& m1) {
        if (m1.size() != m2.size()) {
            throw "bad thing";
        }

        std::vector<T> result(m1.size(), T());
        for (size_t i = 0; i < m1.size(); i++) {
            for (size_t j = 0; j < m1.size(); j++) {
                result[i] += m1[j][i] * m2[j];
            }
        }

        return result;
    }
    std::vector<T>& operator[](const size_t i) {
        return matrix[i];
    }
    std::vector<T> operator[](const size_t i) const {
        return matrix[i];
    }
    friend std::ostream& operator<< (std::ostream& out, const Matrix<T>& matr) {
        for (size_t i = 0; i < matr.size(); i++) {
            for (size_t j = 0; j < matr.size(); j++) {
                out << matr.matrix[i][j] << ' ';
            }
            out << std::endl;
        }
        return out;
    }
};

void Task_1_1() {
    Matrix<double> m(4, {
        1, 2, -1, -7,
        8, 0, -9, -3,
        2, -3, 7, 1,
        2, -5, -6, 8,
    });

    // 0. Исходная матрица
    print("0. Matrix");
    print(m);
    print("Vector b: ");
    print(std::vector<double>{ -23, 39, -7, 30 });
    print();

    // 1. LU - разложение
    print("1. LU - Decomposition");
    auto res = m.LUFactorizing();
    for (auto& mm : res) print(mm);
    print(res[0] * res[1] * res[2]);

    // 2. Решение системы
    print("2. System Solution");
    auto solution = m.Solve({-23, 39, -7, 30});
    print(solution);
    print();

    // 3. Определитель
    print("3. Determinant");
    auto determinant = m.Determinant();
    print(determinant);
    print();

    // 4. Обратная матрица
    print("4. Inverse Matrix");
    auto reverse = m.Reverse();
    print(reverse);
}

int main()
{
    Task_1_1();
    return 0;
}
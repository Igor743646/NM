#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void print() {
    std::cout << std::endl;
}

template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
    for (auto& a : v) {
        out << a << " ";
    }
    return out;
}

template<class T>
void print(T obj) {
    std::cout << obj << std::endl;
}


template<class T>

class Matrix {

    vector<vector<T>> matrix;

public:

    Matrix(size_t _n) {
        matrix = vector<vector<T>>(_n, vector<T>(_n, T()));
    }

    Matrix(size_t _n, vector<T> v) : Matrix(_n) {
        for (size_t i = 0; i < _n; i++) {
            for (size_t j = 0; j < _n; j++) {
                matrix[i][j] = v[i * _n + j];
            }
        }
    }

    size_t inline size() {
        return matrix.size();
    }

    vector<T>& operator[] (size_t i) {
        return matrix[i];
    }

    vector<T> operator[] (size_t i) const {
        return matrix[i];
    }

    Matrix<T> operator* (const Matrix<T>& m) {
        Matrix<T> result(size());

        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                result[i][j] = T();
                for (size_t k = 0; k < size(); k++) {
                    result[i][j] += matrix[i][k] * m[k][j];
                }
            }
        }

        return result;
    }

    friend ostream& operator<<(ostream& out, Matrix<T> m) {
        for (size_t i = 0; i < m.size(); i++) {
            for (size_t j = 0; j < m.size(); j++) {
                out << m[i][j] << " ";
            }
            out << endl;
        }
        return out;
    }
};

template<class T>
T Norm(Matrix<T> A) {
    T sum = T();

    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = i + 1; j < A.size(); j++) {
            sum += A[i][j] * A[i][j];
        }
    }
    return std::sqrt(sum);
}

template<class T>
void JEA(Matrix<T> A, T eps) {

    Matrix<T> Eigenvectors(A.size());
    Matrix<T> U(A.size());
    Matrix<T> U_trans(A.size());

    for (size_t i = 0; i < U.size(); i++) {
        Eigenvectors[i][i] = 1;
        U[i][i] = 1;
        U_trans[i][i] = 1;
    }

    while (Norm(A) > eps) {
        T max = std::abs(A[0][1]);
        size_t l = 0, m = 1;

        // 1. Поиск максимального по модулю недиагонального элемента
        for (size_t i = 0; i < A.size(); i++) {
            for (size_t j = i + 1; j < A.size(); j++) {
                if (std::abs(A[i][j]) > max) {
                    max = std::abs(A[i][j]);
                    l = i;
                    m = j;
                }
            }
        }

        // 2. Вычисление угла поворота
        double phi = 0.5 * (std::atan(2 * A[l][m] / (A[l][l] - A[m][m])));

        // 3. Составление матрицы поворота
        U[l][l] = std::cos(phi);
        U[m][m] = std::cos(phi);
        U[l][m] = -std::sin(phi);
        U[m][l] = std::sin(phi);

        U_trans[l][l] = std::cos(phi);
        U_trans[m][m] = std::cos(phi);
        U_trans[l][m] = std::sin(phi);
        U_trans[m][l] = -std::sin(phi);

        // 4. Поворот
        A = U_trans * A * U;
        Eigenvectors = Eigenvectors * U;

        U[l][l] = 1;
        U[m][m] = 1;
        U[l][m] = 0;
        U[m][l] = 0;

        U_trans[l][l] = 1;
        U_trans[m][m] = 1;
        U_trans[l][m] = 0;
        U_trans[m][l] = 0;
    }

    // Печать результата
    print("Eigenvalues:");
    for (size_t i = 0; i < A.size(); i++) {
        print(A[i][i]);
    }

    print();

    for (size_t j = 0; j < Eigenvectors.size(); j++) {
        for (size_t i = 0; i < Eigenvectors.size(); i++) {
            Eigenvectors[i][j] = Eigenvectors[i][j] / Eigenvectors[Eigenvectors.size() - 1][j];
        }
    }

    print("Eigenvectors:");
    print(Eigenvectors);
}


int main()
{

    Matrix<double> m(3, {
        9, 2, -7,
        2, -4, -1,
        -7, -1, 1
    });

    JEA(m, 0.001);

    return 0;
}
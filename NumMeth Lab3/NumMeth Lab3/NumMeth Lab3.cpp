#include <iostream>
#include <vector>

using namespace std;

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

    void SwapColumns(size_t i, size_t j) {
        for (size_t k = 0; k < size(); k++) {
            std::swap(matrix[k][i], matrix[k][j]);
        }
    }

    void SwapLines(size_t i, size_t j) {
        for (size_t k = 0; k < size(); k++) {
            std::swap(matrix[i][k], matrix[j][k]);
        }
    }

    size_t inline size() {
        return matrix.size();
    }

    vector<T>& operator[] (size_t i) {
        return matrix[i];
    }

    vector<T>& operator[] (size_t i) const {
        return matrix[i];
    }

    vector<T> operator* (vector<T> v) {
        vector<T> result(size(), T());

        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                result[i] += matrix[i][j] * v[j];
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
vector<T> operator+ (vector<T> v1, vector<T> v2) {
    vector<T> result(v1.size(), T());

    for (size_t i = 0; i < v1.size(); i++) {
        result[i] = v1[i] + v2[i];
    }

    return result;
}

template<class T>
vector<T> operator- (vector<T>& v1, vector<T>& v2) {
    vector<T> result(v1.size(), T());

    for (size_t i = 0; i < v1.size(); i++) {
        result[i] = v1[i] - v2[i];
    }

    return result;
}

template<class T>
T Norm(vector<T> v) {
    T max = std::abs(v[0]);

    for (T& e : v) {
        max = std::max(max, std::abs(e));
    }

    return max;
}

template<class T>
T Norm(Matrix<T> A) {
    T max = -1;

    for (size_t i = 0; i < A.size(); i++) {
        T sum = 0;
        for (size_t j = 0; j < A.size(); j++) {
            sum += std::abs(A[i][j]);
        }
        max = std::max(max, std::abs(sum));
    }

    return max;
}

template<class T>
vector<T> Solve(Matrix<T> A, vector<T> b, T eps) {
    vector<T> x(b);
    
    // x = a2 + a1 * x
    Matrix<T> a1(A.size());
    vector<T> a2(A.size(), T());

    // 1. Вычисление a2
    for (size_t i = 0; i < A.size(); i++) {
        a2[i] = b[i] / A[i][i];
    }

    // 2. Вычисление a1
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A.size(); j++) {
            if (i == j) { 
                a1[i][j] = 0; 
            }
            else {
                a1[i][j] = -(A[i][j] / A[i][i]);
            }
        }
    }

    if (Norm(a1) >= 1) {
        print("We have problems, cause a1 matrix's norm >= 1");
        return vector<T>();
    }

    // 3. Итерации
    T a1_norm = Norm(a1);
    vector<T> next_x = a2 + a1 * x;
    T eps_k = (a1_norm) / (1 - a1_norm) * (Norm(next_x - x));

    while (eps_k > eps) {
        x = next_x;
        next_x = a2 + (a1 * next_x);
        eps_k = (a1_norm) / (1 - a1_norm) * (Norm(next_x - x));
    }

    return next_x;
}

int main()
{

    Matrix<double> m1 (4, {
        23, -6, -5, 9,
        8, 22, -2, 5,
        7, -6, 18, -1,
        3, 5, 5, -19,
    });

    auto m2 = vector<double>{232, -82, 202, -57};

    cout << Solve(m1, m2, 0.0001) << endl;

    return 0;
}
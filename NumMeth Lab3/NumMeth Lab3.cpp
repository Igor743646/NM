#include <iostream>
#include <vector>

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

    vector<T> operator[] (size_t i) const {
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

    Matrix<T> operator+ (const Matrix<T>& m) {
        Matrix<T> result(size());

        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                result[i][j] = matrix[i][j] + m[i][j];
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

void DebugSolves(int type) {
    static int count1 = 0;
    static int count2 = 0;

    if (type == 1) {
        count1++;
    }
    else if (type == 2) {
        count2++;
    }
    else if (type == 3) {
        cout << "Yacobi iters: " << count1 << endl;
        cout << "Zeidel iters: " << count2 << endl;
    }
}

template<class T>
vector<T> SolveYacobi(Matrix<T> A, vector<T> b, T eps) {
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
        DebugSolves(1);
        x = next_x;
        next_x = a2 + (a1 * next_x);
        eps_k = (a1_norm) / (1 - a1_norm) * (Norm(next_x - x));
    }

    return next_x;
}

template<class T>
vector<T> SolveZeidel(Matrix<T> A, vector<T> b, T eps) {
    vector<T> x(b);

    // x = a1 * x + a2 * x + b
    Matrix<T> a(A.size());
    Matrix<T> a2(A.size());
    vector<T> new_b(A.size(), T());

    // 1. Вычисление new_b
    for (size_t i = 0; i < A.size(); i++) {
        new_b[i] = b[i] / A[i][i];
    }

    // 2. Вычисление a1 и a2
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A.size(); j++) {
            if (i == j) {
                a[i][j] = 0;
            }
            else {
                a[i][j] = -(A[i][j] / A[i][i]);
            }

            if (i < j) {
                a2[i][j] = -(A[i][j] / A[i][i]);
            }
        }
    }

    if (Norm(a) >= 1) {
        print("We have problems, cause a1 matrix's norm >= 1");
        return vector<T>();
    }

    // 3. Итерации
    T a_norm = Norm(a);
    T a2_norm = Norm(a2);
    vector<T> next_x = new_b + a * x;
    T eps_k = (a2_norm) / (1 - a_norm) * (Norm(next_x - x));

    while (eps_k > eps) {
        DebugSolves(2);
        x = next_x;

        // Зейдельское улучшение
        for (size_t i = 0; i < A.size(); i++) {
            T buf = new_b[i];

            for (size_t j = 0; j < i; j++) {
                buf += a[i][j] * next_x[j];
            }

            for (size_t j = i + 1; j < A.size(); j++) {
                buf += a[i][j] * next_x[j];
            }

            next_x[i] = buf;
        }
        // .....................

        eps_k = (a2_norm) / (1 - a_norm) * (Norm(next_x - x));
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

    cout << SolveYacobi(m1, m2, 0.0001) << endl;
    cout << SolveZeidel(m1, m2, 0.0001) << endl;
    DebugSolves(3);
    return 0;
}
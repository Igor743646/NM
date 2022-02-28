#include <iostream>
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

template<class T>
std::vector<T> Solve(std::vector<std::vector<T>>& abc, std::vector<T> d) {
    size_t dimension = d.size();
    std::vector<T> result(dimension);

    std::vector<T> P(dimension, 0);
    std::vector<T> Q(dimension, 0);
    P[0] = -(abc[2][0] / abc[1][0]);
    Q[0] = (d[0] / abc[1][0]);

    for (size_t i = 1; i < dimension - 1; i++) {
        P[i] = -(abc[2][i] / (abc[1][i] + abc[0][i - 1] * P[i - 1]));
        Q[i] = ((d[i] - abc[0][i - 1] * Q[i - 1]) / (abc[1][i] + abc[0][i - 1] * P[i - 1]));
    }

    result[dimension - 1] = ((d[dimension - 1] - abc[0][dimension - 2] * Q[dimension - 2]) / (abc[1][dimension - 1] + abc[0][dimension - 2] * P[dimension - 2]));

    for (size_t i = 0; i < dimension - 1; i++) {
        size_t k = dimension - 2 - i;
        result[k] = P[k] * result[k + 1] + Q[k];
    }

    return result;
}

void Task_1_2() {
    std::vector<std::vector<double>> abc {
           { -6,   9,  8,   6},
        { 6, 16, -17, 22, -13},
        {-5,  9,  -3, -8}
    };

    std::vector<double> d {-58, 161, -114, -90, -55};

    // 0. Входные данные
    print("0. Entry data");
    for (int i = 0; i < 3; i++) {
        print(std::vector<std::string>{"a:", "b:", "c:"}[i]);
        print(abc[i]);
    }
    print("d:");
    print(d);
    print();

    // 1. Решение
    print("Solution");
    print(Solve(abc, d));
}

int main()
{
    Task_1_2();
    return 0;
}
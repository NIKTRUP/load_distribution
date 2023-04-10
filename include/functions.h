#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <numeric>
#include <cmath>
#include <numbers>
#include <exception>
#include <string>
#include <cassert>
#include <utility>
#include <fstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include "include/write.h"


using ldouble = long double;
using namespace std::literals;

// Эта ошибка должна выбрасываться при ошибке в входных параметрах
class ArgumentError : public std::invalid_argument {
public:
    using invalid_argument::invalid_argument;
};

// Эта ошибка должна выбрасываться при ошибке в вычислениях
class CalculateError : public std::runtime_error {
public:
    using runtime_error::runtime_error;
};

// Эта ошибка должна выбрасываться при отсутствии решения для СЛАУ
class NoSolution : public std::runtime_error {
public:
    using runtime_error::runtime_error;
};

// проверяет на равенство два числа ldouble с точностью до 15 знаков после запятой
bool is_equal(ldouble x, ldouble y);

// Проверяет принадлежность точки x к отрезку
bool in_segment(ldouble x, ldouble begin, ldouble end);

// Перегрузка
bool in_segment(ldouble x, std::pair<ldouble, ldouble> seg);

// создаёт равномерно распределенную последовательность в указанном интервале.
template <typename T>
std::vector<T> linspace(T start, T end, ldouble num){
    std::vector<T> linspaced;
    if (0 != num){
        if (1 == num) {
            linspaced.push_back(start);
        }
        else if(2 == num){
            linspaced.push_back(start);
            linspaced.push_back(end);
        }
        else{
            ldouble delta = (end - start) / (num - 1);
            for (auto i = 0; i < (num - 1); ++i){
                linspaced.push_back(start + delta * i);
            }
            linspaced.push_back(end);
        }
    }
    return linspaced;
}


// Прямой ход метода Гаусса
template <typename Matrix, typename  Vector>
void Triangulation(Matrix& matrix, Vector& b){
    ldouble max, s;
    size_t size = matrix.size();
    for (size_t k = 0; k < size - 1; ++k) {
        max = std::abs(matrix[k][k]);
        size_t n = k;
        for (size_t i = k + 1; i < size ; ++i) {
            s = std::abs(matrix[i][k]);
            if (s > max) {
                max = s;
                n = i;
            }
        }
        if (n > k) {
            for (size_t j = k; j< size ; ++j) {
                std::swap(matrix[k][j], matrix[n][j]);
            }
            std::swap(b[k], b[n]);
        }
        for (size_t i= k + 1; i < size ; ++i) {
            auto mul = matrix[i][k]/matrix[k][k];
            for (size_t j=k+1; j < size ; ++j) {
                matrix[i][j] = matrix[i][j] - matrix[k][j]*mul;
            }
            b[i]=b[i]-b[k]*mul;
        }
    }
}

// Решение СЛАУ Методом Гаусса
template <typename Matrix, typename  Vector>
void SolveGauss(Matrix& matrix, Vector& b) {
    size_t size = matrix.size();
    Triangulation(matrix, b);

    b[size - 1] = b[size - 1]/matrix[size - 1][size - 1];
    for (size_t i = 0; i < size - 1; ++i) {
        size_t k = size-(i+1)-1;
        for (size_t j=k+1; j < size; ++j){
            b[k] -= matrix[k][j]*b[j];
        }
        b[k] /= matrix[k][k];
    }
}
#endif // FUNCTIONS_H

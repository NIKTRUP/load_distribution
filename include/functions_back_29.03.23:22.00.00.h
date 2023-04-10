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

// Эта ошибка должна выбрасываться при ошибке в вычислениях
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

// если точка x находится в отрезке от begin до end, то true, иначе false
bool in_segment(ldouble x, ldouble begin, ldouble end);

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

// T0 - натяжение ненагруженной струны
// Функция Грина
ldouble G(ldouble x, ldouble xi, ldouble T0 ,ldouble begin, ldouble end);

// Функция H
ldouble H(ldouble xi, ldouble x, std::vector<ldouble> t ,ldouble T0 ,ldouble begin, ldouble end);

// Функция H Янвная
ldouble H(ldouble xi, ldouble x, ldouble T0 , ldouble begin, ldouble end);

template <typename T>
struct Points{
  const std::vector<T>& x;
  const std::vector<T>& y;
};

// Вычисляет отклонение от положения равновесия
ldouble СalculateDeviationFromEquilibrium(ldouble x, Points<ldouble> p,
                                         ldouble T0 , ldouble begin, ldouble end);

// Вычисляет сумму квадрата ошибки между заданным отклонением и вычисленным
ldouble CalculateLoss(std::vector<ldouble> x, std::vector<ldouble> u, std::vector<ldouble> u0);

// Выводит СЛАУ
template <typename Matrix, typename Vector>
void OutSLAE(const Matrix& matrix, const Vector& free_term_column){
    size_t n = free_term_column.size();
    for (size_t i = 0; i < n; ++i){
      for (size_t j = 0; j < n; ++j){
          std::cout << matrix[i][j] << "*x" << j;
        if (j < n - 1){
            std::cout << " + ";
        }
      }
      std::cout << " = " << free_term_column[i] << '\n';
    }
    return;
}

ldouble GetStep(const std::vector<ldouble>& x, size_t i);

// Заполняет СЛАУ
template <typename Matrix, typename Vector>
void FillSLAE(Matrix& matrix, Vector& b, const std::vector<ldouble>& x,
              const std::vector<ldouble>& u0,ldouble alpha, ldouble P, ldouble T0,
              ldouble begin, ldouble end){
    auto n = matrix.size();
    ldouble h = 0.0;
    for(size_t i = 0; i < n; ++i){

        b[i] = (i == n-1) ? P : СalculateDeviationFromEquilibrium(x[i], {x, u0}, T0, begin, end); // заполняю правую часть
        matrix[i][n-1] = (i == n-1) ? 0 : 0.5; // заполняю последний столбец матрицы

         //h = x[1] - x[0];
        for(size_t j = 0; j < n - 1; ++j){

            h = GetStep(x, j);

            matrix[i][j] = 0;
            if(i == j){
                matrix[i][j] += alpha;
            }

            if(i == n - 1){ // заполняем последнюю строку
                matrix[i][j] = h/2.0;
            }else { // заполняем от первого до предпоследнего стобеца (кроме последней строки)
                matrix[i][j] += H(x[j], x[i], T0, begin, end)*h/2.0;
            }

        }
    }
}

//if(i == n - 1){ // заполняем последнюю строку
//    matrix[i][j] = (j == 0 || j == n - 2) ? h/2 : h;
//}else { // заполняем от первого до предпоследнего стобеца (кроме последней строки)
//    matrix[i][j] += H(x[j], x[i], T0, begin, end)*h/2.0;
//}

//if(i == n - 1){ // заполняем последнюю строку
//    matrix[i][j] = (j == 0 || j == n - 2) ? h/2 : h;
//}else if(i == 0 || i == n - 2){ // заполняем от первого до предпоследнего стобеца (кроме последней строки)
//    matrix[i][j] += H(x[j], x[i], T0, begin, end)*h/2.0;
//}else{
//    matrix[i][j] += H(x[j], x[i], T0, begin, end)*h;
//}

// Приводит матрицу к верхней треугольной
template <typename Matrix, typename Vector>
void Triangulation(Matrix& matrix, Vector& free_term_column, size_t n){
    ldouble eps = 1e-33;
    ldouble max;

    size_t k = 0, index;
    while (k < n){
      max = std::abs(matrix[k][k]);
      index = k;
      for (size_t i = k + 1; i < n; i++){
        if (std::abs(matrix[i][k]) > max){
          max = std::abs(matrix[i][k]);
          index = i;
        }
      }
      // Перестановка строк
      if (max <= eps) {
          WriteSLAE(matrix, free_term_column, "slae_no_solution");
          throw NoSolution("Решение получить невозможно из-за нулевого столбца " + std::to_string(index) + " матрицы");
      }
      for (size_t j = 0; j < n; j++){
        ldouble temp = matrix[k][j];
        matrix[k][j] = matrix[index][j];
        matrix[index][j] = temp;
      }
      ldouble temp = free_term_column[k];
      free_term_column[k] = free_term_column[index];
      free_term_column[index] = temp;

      for (size_t i = k; i < n; i++) {
        ldouble temp = matrix[i][k];
        if (std::abs(temp) <= eps) { continue; }
        for (size_t j = 0; j < n; j++){
           matrix[i][j] = matrix[i][j] / temp;
        }
        free_term_column[i] = free_term_column[i] / temp;
        if (i == k)  { continue; }
        for (size_t j = 0; j < n; j++){
            matrix[i][j] = matrix[i][j] - matrix[k][j];
        }
        free_term_column[i] = free_term_column[i] - free_term_column[k];
      }
      k++;
    }
}

// Решает СЛАУ методом Гаусса
template <typename Matrix, typename Vector>
Vector SolveGaussMethod(Matrix& matrix, Vector& free_term_column) {
    size_t n = free_term_column.size();
    Vector solution(n);

    Triangulation(matrix, free_term_column, n);

    for (int j = n - 1; j >= 0; --j){
       solution[j] = free_term_column[j];
       for (int i = 0; i < j; i++){
           free_term_column[i] -= matrix[i][j] * solution[j];
       }
    }
    return solution;
}


// СЛАУ, метод Гаусса
template <typename Matrix, typename  Vector>
void gauss(Matrix& a, Vector& b) {
  ldouble z, s;
  int min = 0; int max=a.size() - 1;
  for (int k=min; k<=max-1; ++k) {
    z=fabs(a[k][k]);
    int n=k;
    for (int i=k+1; i<=max; ++i) {
      s=fabs(a[i][k]);
      if (s>z) {z=s; n=i;};
      };
    if (n>k) {
      for (int j=k; j<=max; ++j) {
        z=a[k][j]; a[k][j]=a[n][j]; a[n][j]=z;
        };
      z=b[k]; b[k]=b[n]; b[n]=z;
      };
    for (int i=k+1; i<=max; ++i) {
      z=a[i][k]/a[k][k];
      for (int j=k+1; j<=max; ++j) a[i][j]=a[i][j]-a[k][j]*z;
      b[i]=b[i]-b[k]*z;
      };
    };
  b[max]=b[max]/a[max][max];
  for (int i=min; i<=max-1; ++i) {
    int k=max-(i-min+1);
    for (int j=k+1; j<=max; ++j) b[k]=b[k]-a[k][j]*b[j];
    b[k]=b[k]/a[k][k];
    };
  return;
  };

#endif // FUNCTIONS_H

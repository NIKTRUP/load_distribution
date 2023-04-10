#ifndef STRING_LOAD_DISTRIBUTION_H
#define STRING_LOAD_DISTRIBUTION_H
#include "functions.h"

namespace str_dstrb {
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

             h = x[1] - x[0];
            for(size_t j = 0; j < n - 1; ++j){

                //h = GetStep(x, j);

                matrix[i][j] = 0;
                if(i == j){
                    matrix[i][j] += alpha;
                }

                if(i == n - 1){ // заполняем последнюю строку
                    matrix[i][j] = h; // h/2.0
                }else { // заполняем от первого до предпоследнего стобеца (кроме последней строки)
                    matrix[i][j] += H(x[j], x[i], T0, begin, end)*h; // h/2.0
                }

            }
        }
    }
}

#endif // STRING_LOAD_DISTRIBUTION_H

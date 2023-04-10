#ifndef MEMBRANE_LOAD_DISTRIBUTION_H
#define MEMBRANE_LOAD_DISTRIBUTION_H
#include "include/functions.h"
using ldouble = long double;

namespace membrane_dstrb {

    template <typename T>
    struct Coordinates{
        const std::vector<T>& x;
        const std::vector<T>& y;
    };

    template <typename T>
    struct Points{
      const std::vector<T>& x;
      const std::vector<T>& y;
      const std::vector<std::vector<T>>& z;
    };

    ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y, ldouble end_x, ldouble end_y, size_t M, size_t N);

    ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y, ldouble end_x, ldouble end_y, ldouble eps_x, ldouble eps_y);

    ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y, ldouble end_x, ldouble end_y, ldouble eps);

    namespace detail {
        ldouble СalculateDeviationFromEquilibrium(ldouble x, ldouble y, Points<ldouble> p,
                                                  size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                                                  ldouble end_x, ldouble end_y, bool set_limit);

        std::vector<std::vector<ldouble>> CalculateU(const Points<ldouble> p, size_t M, size_t N, ldouble eps_x,
                                                     ldouble eps_y, size_t end_x, size_t end_y, bool set_limit);

        ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
                  ldouble x, ldouble y, size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, bool set_limit);

        ldouble H(ldouble xi, ldouble eta,
                  ldouble x, ldouble y, size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, bool set_limit);

        ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y,
                  ldouble end_x, ldouble end_y, size_t M, size_t N,
                  ldouble eps_x, ldouble eps_y, bool set_limit);
    }

    ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
                      ldouble x, ldouble y, size_t M, size_t N,
                      ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps_x, ldouble eps_y,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    ldouble H(ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps_x, ldouble eps_y,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    ldouble H(ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    ldouble H(ldouble xi, ldouble eta,
              ldouble x, ldouble y, size_t M, size_t N,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    // Вычисляет отклонение от положения равновесия
    ldouble СalculateDeviationFromEquilibrium(ldouble x, ldouble y, Points<ldouble> p,
                                              size_t M, size_t N, ldouble end_x, ldouble end_y);

    ldouble СalculateDeviationFromEquilibrium(ldouble x, ldouble y, Points<ldouble> p,
                                              ldouble eps_x, ldouble eps_y, ldouble end_x, ldouble end_y);

    // Вычисляет сумму квадрата ошибки между заданным отклонением и вычисленным
    ldouble CalculateLoss(const Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u,
                          const std::vector<std::vector<ldouble>>& u0);

    std::vector<std::vector<ldouble>> CalculateU(const Points<ldouble> p, size_t M, size_t N,
                                                 size_t end_x, size_t end_y);

    std::vector<std::vector<ldouble>> CalculateU(const Points<ldouble> p, ldouble eps_x,
                                                 ldouble eps_y, size_t end_x, size_t end_y);

    // Заполняет СЛАУ
    template <typename Matrix, typename Vector>
    void FillSLAE(Matrix& matrix, Vector& b, const Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u0, ldouble alpha, ldouble P, size_t size_x, size_t size_y, size_t M, size_t N,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y){

        ldouble hx = p.x[1] - p.x[0],
                hy = p.y[1] - p.y[0];


        auto H_ = [&](ldouble xi, ldouble eta, ldouble x, ldouble y){
            return H(p, xi, eta, x, y, M, N, begin_x, end_x, begin_y, end_y);
        };


        auto fill_block = [&](size_t k){
            for(size_t i = 0; i < size_y; ++i){
                for(size_t j = 0; j < size_x; ++j){
                    matrix[i + k*size_y][j + k*size_x] += H_(p.x[j], p.y[i], p.x[k], p.y[i])*hx*hy;
                }
                b[i + k*size_y] = СalculateDeviationFromEquilibrium(p.x[k], p.y[i], {p.x, p.y, u0}, M, N, end_x, end_y);
            }
        };

        // Заполняет блоки и правую часть
        for(size_t k = 0; k < size_y; ++k){
            fill_block(k);
        }
        b[size_y*size_x] = P;

        // Заполняет последнюю строку и столбец. Добавляет параметр регуляризации к диагонали;
        for(size_t i = 0; i < size_x*size_y; ++i){
            matrix[size_x*size_y][i] = hx*hy;
            matrix[i][size_x*size_y] = 0.5;
            matrix[i][i] += alpha;
        }
    }

    std::vector<std::vector<ldouble>> FillP(std::vector<ldouble>& b, size_t size_x, size_t size_y);

    std::vector<std::vector<ldouble>> AddZerosBorders(std::vector<std::vector<ldouble>>& u);
}


#endif // MEMBRANE_LOAD_DISTRIBUTION_H

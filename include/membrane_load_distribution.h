#ifndef MEMBRANE_LOAD_DISTRIBUTION_H
#define MEMBRANE_LOAD_DISTRIBUTION_H
#include "functions.h"
#include <omp.h>
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


    ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
                      ldouble x, ldouble y, size_t M, size_t N,
                      ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    ldouble H(Coordinates<ldouble> p, ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps_x, ldouble eps_y,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y);

    ldouble H(Coordinates<ldouble> p, ldouble xi, ldouble eta,
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

    namespace detail {
        ldouble СalculateDeviationFromEquilibrium(ldouble x, ldouble y, Points<ldouble> p,
                                                  size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                                                  ldouble end_x, ldouble end_y, bool set_limit);

        std::vector<std::vector<ldouble>> CalculateU(Points<ldouble> p, size_t M, size_t N, ldouble eps_x,
                                                     ldouble eps_y, ldouble end_x, ldouble end_y, bool set_limit);

        ldouble H(Coordinates<ldouble> p, ldouble xi, ldouble eta,
                  ldouble x, ldouble y, size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, bool set_limit);

        ldouble H(ldouble xi, ldouble eta,
                  ldouble x, ldouble y, size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, bool set_limit);

        ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y,
                  ldouble end_x, ldouble end_y, size_t M, size_t N,
                  ldouble eps_x, ldouble eps_y, bool set_limit);

        template <typename Matrix, typename Vector>
        void FillSLAE(Matrix& matrix, Vector& b, const Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u0,
                      ldouble alpha, ldouble P, size_t size_x, size_t size_y, size_t M, size_t N, ldouble eps_x, ldouble eps_y, bool set_limit,
                      ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, size_t num_threads) {

            ldouble hx = p.x[1] - p.x[0],
                    hy = p.y[1] - p.y[0];


            if (size_x < 5 && size_y < 5) {
                num_threads = 1;
            }

            auto H_ = [&](ldouble xi, ldouble eta, ldouble x, ldouble y) {
                if (set_limit) {
                    return membrane_dstrb::H(xi, eta, x, y, M, N, begin_x, end_x, begin_y, end_y);
                } else {
                    return membrane_dstrb::H(xi, eta, x, y, eps_x, eps_y, begin_x, end_x, begin_y, end_y);
                }
            };

            auto CalculateDeviation = [&](ldouble xi, ldouble yk) {
                if (set_limit) {
                    return membrane_dstrb::СalculateDeviationFromEquilibrium(xi, yk, {p.x, p.y, u0}, M, N, end_x,
                                                                             end_y);
                } else {
                    return membrane_dstrb::СalculateDeviationFromEquilibrium(xi, yk, {p.x, p.y, u0}, eps_x, eps_y,
                                                                             end_x, end_y);
                }
            };

            auto fill_block = [&](size_t m, size_t k) {
                #pragma omp parallel for
                for (size_t i = 0; i < size_x; ++i) {
                    #pragma omp parallel for
                    for (size_t j = 0; j < size_x; ++j) {
                        matrix[i + m * size_x][j + k * size_x] = H_(p.x[j], p.y[k], p.x[i], p.y[m]) * hx * hy;
                    }
                }
            };

            omp_set_num_threads(num_threads);
            #pragma omp parallel for
            for (size_t m = 0; m < size_y; ++m) {
                #pragma omp parallel for
                for (size_t k = 0; k < size_y; ++k) {
                    fill_block(m, k);
                }
                #pragma omp parallel for
                for(size_t i = 0; i < size_x; ++i){
                    b[i + m * size_x] = CalculateDeviation(p.x[i], p.y[m]);
                }
            }
            b[size_y * size_x] = P;

            // Заполняет последнюю строку и столбец. Добавляет параметр регуляризации к диагонали;
            for (size_t i = 0; i < size_x * size_y; ++i) {
                matrix[size_x * size_y][i] = hx * hy;
                matrix[i][size_x * size_y] = 0.5;
                matrix[i][i] += alpha;
            }
        }
    }

    // Вычисляет сумму квадрата ошибки между заданным отклонением и вычисленным
    ldouble CalculateLoss(Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u,
                          const std::vector<std::vector<ldouble>>& u0);

    std::vector<std::vector<ldouble>> CalculateU(Points<ldouble> p, size_t M, size_t N,
                                                 ldouble end_x, ldouble end_y);

    std::vector<std::vector<ldouble>> CalculateU(Points<ldouble> p, ldouble eps_x,
                                                 ldouble eps_y, ldouble end_x, ldouble end_y);

    std::vector<std::vector<ldouble>> CalculateU(Points<ldouble> p, ldouble eps,
                                                 ldouble end_x, ldouble end_y);

    // Заполняет СЛАУ
    template <typename Matrix, typename Vector>
    void FillSLAE(Matrix& matrix, Vector& b, const Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u0,
                  ldouble alpha, ldouble P, size_t size_x, size_t size_y, size_t M, size_t N,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, size_t num_threads = 6){
        return detail::FillSLAE(matrix, b, p, u0, alpha, P, size_x, size_y, M, N, 0.0, 0.0, true, begin_x, end_x, begin_y, end_y, num_threads);
    }

    template <typename Matrix, typename Vector>
    void FillSLAE(Matrix& matrix, Vector& b, const Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u0,
                  ldouble alpha, ldouble P, size_t size_x, size_t size_y, ldouble eps_x, ldouble eps_y,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, size_t num_threads = 6){
        return detail::FillSLAE(matrix, b, p, u0, alpha, P, size_x, size_y, 0, 0, eps_x, eps_y, false, begin_x, end_x, begin_y, end_y, num_threads);
    }

    template <typename Matrix, typename Vector>
    void FillSLAE(Matrix& matrix, Vector& b, const Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u0,
                  ldouble alpha, ldouble P, size_t size_x, size_t size_y, ldouble eps,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, size_t num_threads = 6){
        return detail::FillSLAE(matrix, b, p, u0, alpha, P, size_x, size_y, 0, 0, eps, eps, false, begin_x, end_x, begin_y, end_y, num_threads);
    }

    std::vector<std::vector<ldouble>> FillP(std::vector<ldouble>& b, size_t size_x, size_t size_y);

    std::vector<std::vector<ldouble>> AddZerosBorders(std::vector<std::vector<ldouble>>& u);
}

#endif // MEMBRANE_LOAD_DISTRIBUTION_H

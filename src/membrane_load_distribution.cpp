#include "include/membrane_load_distribution.h"

namespace membrane_dstrb {

    namespace detail {
        ldouble СalculateDeviationFromEquilibrium(ldouble x, ldouble y, Points<ldouble> p,
                                                  size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                                                  ldouble end_x, ldouble end_y, bool set_limit){
            const size_t & a = end_x,  // по строкам: x -> b
                         & b = end_y; // по столбцам: y -> a

            ldouble pi_4 = M_PI*M_PI*M_PI*M_PI,
                    mul = 4.0/(pi_4*b*a),
                    hy = p.y[1] - p.y[0],
                    hx = p.x[1] - p.x[0];
            ldouble u = 0.0;

            auto G_ = [&](ldouble xj, ldouble yi){
                if(set_limit){
                    return membrane_dstrb::G(xj, yi, x, y, a, b, M, N);
                }else{
                    return membrane_dstrb::G(xj, yi, x, y, a, b, eps_x, eps_y);
                }
            };

            for(size_t i = 0; i < p.y.size(); ++i){
                for(size_t j = 0; j < p.x.size(); ++j){
                    u += G_(p.x[j], p.y[i])*p.z[j][i]*hy*hx;
                }
            }
            u *= mul;
            return  u;
        }

        std::vector<std::vector<ldouble>> CalculateU(const Points<ldouble> p, size_t M, size_t N, ldouble eps_x,
                                                     ldouble eps_y, size_t end_x, size_t end_y, bool set_limit){

            size_t size_x = p.x.size(), size_y = p.y.size();
            std::vector<std::vector<ldouble>> u(size_y + 2, std::vector<ldouble>(size_x + 2));
            ldouble hx = p.x[1] - p.x[0],
                    hy = p.y[1] - p.y[0];
            u[0] = std::vector<ldouble>(size_x + 2, 0);
            for(size_t i = 1 ; i < size_y + 1; ++i){
                u[i][0] = 0;
                for(size_t j = 1; j < size_x + 1 ; ++j){
                    u[i][j] = membrane_dstrb::detail::СalculateDeviationFromEquilibrium(p.x[j - 1], p.y[i - 1], p, M, N,
                                                                                        eps_x, eps_y, end_x, end_y, set_limit)*hx*hy;
                }
                u[i].back() = 0;
            }
            u.back() = std::vector<ldouble>(size_x + 2, 0);
            return u;
        }

        ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
                  ldouble x, ldouble y, size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, bool set_limit){

            auto is_inside_segment = [](ldouble x, const std::string& name,  ldouble begin, ldouble end){
                if(!in_segment(x, begin, end)){
                    throw ArgumentError(name + " = " + std::to_string(x) + " находится не в интервале от " + std::to_string(begin) + " до " + std::to_string(end));
                }
            };

            is_inside_segment(xi, "xi", begin_x, end_x);
            is_inside_segment(x, "x", begin_x, end_x);

            is_inside_segment(eta, "eta", begin_y, end_y);
            is_inside_segment(y, "y", begin_y, end_y);

            auto G_ = [&](ldouble xi,ldouble eta,ldouble x,ldouble y){
                if(set_limit){
                    return membrane_dstrb::G(xi, eta, x, y, end_x, end_y, M, N);
                }else{
                    return membrane_dstrb::G(xi, eta, x, y, end_x, end_y, eps_x, eps_y);
                }
            };

            ldouble sum = 0,
                    hx = p.x[1] - p.x[0],
                    hy = p.y[1] - p.y[0];
            for(size_t i = 0 ; i < p.y.size(); ++i){
                for(size_t j = 0; j < p.x.size(); ++j){
                    sum += G_(x, y, p.x[j], p.y[i])*G_(xi, eta, p.x[j], p.y[i])*hx*hy;
                }
            }
            return sum;
        }

        ldouble H(ldouble xi, ldouble eta,
                  ldouble x, ldouble y, size_t M, size_t N, ldouble eps_x, ldouble eps_y,
                  ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y, bool set_limit){

            auto is_inside_segment = [](ldouble x, const std::string& name,  ldouble begin, ldouble end){
                if(!in_segment(x, begin, end)){
                    throw ArgumentError(name + " = " + std::to_string(x) + " находится не в интервале от " + std::to_string(begin) + " до " + std::to_string(end));
                }
            };

            is_inside_segment(xi, "xi", begin_x, end_x);
            is_inside_segment(x, "x", begin_x, end_x);

            is_inside_segment(eta, "eta", begin_y, end_y);
            is_inside_segment(y, "y", begin_y, end_y);

            ldouble &a = end_x,
                    &b = end_y,
                    aa = a*a,
                    bb = b*b;

            ldouble H = 0.0, h_n = 0.0, h_m = 0.0;
            auto check_m = [&](size_t m){
                if(set_limit){
                    return m <= M;
                }else{
                    return m == 1 || (std::abs(h_m - H) > eps_y);
                }
            };

            auto check_n = [&](size_t n){
                if(set_limit){
                    return n <= N;
                }else{
                    return n == 1 || (std::abs(h_n - H) > eps_x);
                }
            };

            for(size_t m = 1 ; check_m(m); ++m){
                h_m = H;
                ldouble m_pi_a = m*M_PI/a,
                        sin_m_pi_a_x = std::sin(m_pi_a*x),
                        sin_m_pi_a_xi = std::sin(m_pi_a*xi);
                for(size_t n = 1; check_n(n); ++n){
                    h_n = H;
                    ldouble divider = ((m*m)/(aa) + (n*n)/(bb)),
                    n_pi_b = n*M_PI/b;
                    H += (sin_m_pi_a_x*std::sin(n_pi_b*y)*sin_m_pi_a_xi*std::sin(n_pi_b*eta))/(divider*divider*divider*divider);
                }
            }
            H *= (a*b)/4.0;
            return H;
        }


        ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y,
                  ldouble end_x, ldouble end_y, size_t M, size_t N,
                  ldouble eps_x, ldouble eps_y, bool set_limit){
            if(xi > end_x || y > end_y){
                throw ArgumentError("Выход за границы");
            }

            ldouble g = 0.0, h_m = 0.0, h_n = 0.0,
                    &a = end_x,
                    &b = end_y,
                    aa = a*a,
                    bb = b*b;

            auto check_m = [&](size_t m){
                if(set_limit){
                    return m <= M;
                }else{
                    return m == 1 || std::abs(h_m - g) > eps_y;
                }
            };

            auto check_n = [&](size_t n){
                if(set_limit){
                    return n <= N;
                }else{
                    return n == 1 || std::abs(h_n - g) > eps_x;
                }
            };


            for(size_t m = 1 ; check_m(m); ++m){
                h_m = g;
                ldouble m_pi_a = m*M_PI/a,
                        sin_m_x = std::sin(m_pi_a*x),
                        sin_m_xi = std::sin(m_pi_a*xi);
                for(size_t n = 1; check_n(n); ++n){
                    h_n = g;
                    ldouble divider = ((m*m)/(aa) + (n*n)/(bb));
                    g += ( sin_m_x*std::sin(n*M_PI*y/b)*sin_m_xi*std::sin(n*M_PI*eta/b) )/(divider*divider);
                }
            }
            return g;
        }
    }

    ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y, ldouble end_x, ldouble end_y, size_t M, size_t N){
        return detail::G(xi, eta, x, y, end_x, end_y, M, N, 0.0, 0.0, true);
    }

    ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y, ldouble end_x, ldouble end_y, ldouble eps_x, ldouble eps_y){
        return detail::G(xi, eta, x, y, end_x, end_y, 0, 0, eps_x, eps_y, false);
    }

    ldouble G(ldouble xi, ldouble eta, ldouble x, ldouble y, ldouble end_x, ldouble end_y, ldouble eps){
        return G(xi, eta, x, y, end_x, end_y, eps, eps);
    }

    ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
              ldouble x, ldouble y, size_t M, size_t N,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y){
        return detail::H(p, xi, eta, x, y, M, N, 0.0, 0.0, begin_x, end_x, begin_y, end_y, true);
    }

    ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps_x, ldouble eps_y,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y){
        return detail::H(p, xi, eta, x, y, 0, 0, eps_x, eps_y, begin_x, end_x, begin_y, end_y, false);
    }

    ldouble H(const Coordinates<ldouble> p, ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y){
        return detail::H(p, xi, eta, x, y, 0, 0, eps, eps, begin_x, end_x, begin_y, end_y, false);
    }

    ldouble H(ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps_x, ldouble eps_y,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y){
        return detail::H(xi, eta, x, y, 0, 0, eps_x, eps_y, begin_x,  end_x,  begin_y,  end_y, false);
    }

    ldouble H(ldouble xi, ldouble eta,
              ldouble x, ldouble y, ldouble eps,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y){
        return detail::H(xi, eta, x, y, 0, 0, eps, eps, begin_x,  end_x,  begin_y,  end_y, false);
    }

    ldouble H(ldouble xi, ldouble eta,
              ldouble x, ldouble y, size_t M, size_t N,
              ldouble begin_x, ldouble end_x, ldouble begin_y, ldouble end_y){
        return detail::H(xi, eta, x, y, M, N, 0.0, 0.0, begin_x,  end_x,  begin_y,  end_y, true);
    }


    // Вычисляет отклонение от положения равновесия
    ldouble СalculateDeviationFromEquilibrium(ldouble x, ldouble y, Points<ldouble> p,
                                              size_t M, size_t N, ldouble end_x, ldouble end_y){
        return detail::СalculateDeviationFromEquilibrium(x, y, p, M, N, 0.0, 0.0, end_x, end_y, true);
    }


    // Вычисляет отклонение от положения равновесия
    ldouble СalculateDeviationFromEquilibrium(ldouble x, ldouble y, Points<ldouble> p,
                                              ldouble eps_x, ldouble eps_y, ldouble end_x, ldouble end_y){
        return detail::СalculateDeviationFromEquilibrium(x, y, p, 0.0, 0.0, eps_x, eps_y, end_x, end_y, false);
    }

    // Вычисляет сумму квадрата ошибки между заданным отклонением и вычисленным
    ldouble CalculateLoss(const Coordinates<ldouble> p, const std::vector<std::vector<ldouble>>& u,
                          const std::vector<std::vector<ldouble>>& u0){
        ldouble sum = 0.0,
                loss = 0.0,
                hx = p.x[1] - p.x[0],
                hy = p.y[1] - p.y[0];
        for(size_t i = 0 ; i < p.y.size(); ++i){
            for(size_t j = 0; j < p.x.size() ; ++j){
                loss = (u[i][j] - u0[i][j]);
                sum += (loss*loss)*hx*hy;
            }
        }
        return sum;
    }

    std::vector<std::vector<ldouble>> CalculateU(const Points<ldouble> p, size_t M,
                                                 size_t N, size_t end_x, size_t end_y){
        return detail::CalculateU(p, M, N, 0.0, 0.0, end_x, end_y, true);
    }

    std::vector<std::vector<ldouble>> CalculateU(const Points<ldouble> p, ldouble eps_x,
                                                 ldouble eps_y, size_t end_x, size_t end_y){
        return detail::CalculateU(p, 0, 0, eps_x, eps_y, end_x, end_y, false);
    }

    std::vector<std::vector<ldouble>> FillP(std::vector<ldouble>& b, size_t size_x, size_t size_y){
        std::vector<std::vector<ldouble>> p(size_y);
        for(size_t i = 0; i < p.size(); ++i){
            p[i].insert(
                  p[i].begin(),
                  std::make_move_iterator(b.begin() + i*size_x),
                  std::make_move_iterator(b.begin() + (i+1)*size_x)
                );
        }
        return p;
    }

    std::vector<std::vector<ldouble>> AddZerosBorders(std::vector<std::vector<ldouble>>& u){

        size_t size_y = u.size();
        size_t size_x = u[0].size();

        std::vector<std::vector<ldouble>> out(size_y+2);
        out[0] = std::vector<ldouble>(size_x + 2, 0);
        for(size_t i = 1; i < out.size() - 1; ++i){
            out[i].push_back(0);
            out[i].insert(
                  out[i].begin() + 1,
                  std::make_move_iterator(u[i - 1].begin()),
                  std::make_move_iterator(u[i - 1].end())
                );
            out[i].push_back(0);
        }
        out.back() = std::vector<ldouble>(size_x + 2, 0);
        return out;
    }
}

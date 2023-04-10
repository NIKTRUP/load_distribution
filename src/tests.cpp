#include "include/tests.h"
#include <string_view>
#include <iostream>
namespace tests{
    namespace detail {
        void GenerateLinspace(ldouble begin, ldouble end, size_t size_dots){
            auto x = linspace(begin, end, size_dots);
            ldouble h = x[1] - x[0];
            for(size_t i = 0 ; i < x.size(); ++i){
                bool check_segment = in_segment(x[i], begin, end);
                if(!check_segment){
                    std::cout << "x["<< i << "] = " <<x[i] << '\n';
                }
                ASSERT(check_segment);
                if(i > 0){
                    bool check_h = is_equal(x[i] - x[i - 1], h);
                    if(!check_h){
                        std::cout << "h = " << h << '\n';
                        ldouble real_h = x[i] - x[i-1];
                        std::cout << "x["<< i << "] -" << "x["<< i - 1 << "] = " << real_h << '\n';
                        std::cout << "h - real_h = " << h - real_h << '\n';
                    }
                    ASSERT(is_equal(x[i] - x[i - 1], h));
                }
            }
        }

        std::vector<ldouble> SplitIntoNums(std::string_view str, std::string delimiter) {
            std::vector<ldouble> result;

            while (!str.empty()) {
                size_t space = str.find(delimiter);
                auto str_view = space == str.npos ? str.substr(0) : str.substr(0, space);
                result.push_back(std::stold(std::string(str_view)));
                str.remove_prefix(str_view.size());
                str.remove_prefix(std::min(str.find_first_not_of(delimiter), str.size()));
            }
            return result;
        }

        std::tuple<std::vector<std::vector<ldouble>>, std::vector<ldouble>> GenerateSLAE(size_t size){
            std::vector<std::vector<ldouble>> matrix(size, std::vector<ldouble>(size));
            std::vector<ldouble> b(size);
            for(size_t i = 0 ; i < size; ++i){
                b[i] = std::rand();
                for(size_t j = 0 ; j < size; ++j){
                    matrix[i][j] = std::rand();
                }
            }
            return {matrix, b};
        }
    }

    void TestLinspace(){
        {
            detail::GenerateLinspace(0, 1, 1000);
        }

        {
            detail::GenerateLinspace(0, 1, 100000);
        }

        {
            detail::GenerateLinspace(0, M_PI, 1000);
        }

        {
            detail::GenerateLinspace(0, M_PI, 100000);
        }

        {
            detail::GenerateLinspace(-M_PI, M_PI, 100000);
        }
    }

    void TestIsEqual(){
        ASSERT(is_equal(0.0, 0.0));
        ASSERT(is_equal(M_PI, M_PI));
    }

    void TestIsSegment(){
        ASSERT(in_segment(0.0, 0.0, 1.0));
        ASSERT(in_segment(1.0, 0.0, 1.0));
        ASSERT(in_segment(M_PI, 0.0, M_PI));
        ASSERT(in_segment(3.1415, 0.0, M_PI));

        ASSERT(!in_segment(-1, 0.0, 1.0));
        ASSERT(!in_segment(2, 0.0, 1.0));
    }


    namespace str_dstrb {

        using namespace ::str_dstrb;

        void TestGrinFunction(){
            ldouble eps = 1e-12;
            ldouble begin = 0, end = 1, T0 = 1, x = 0.4, xi = 0.5;
            ASSERT(is_equal(0.2 , G(x, xi, T0, begin, end)));

            begin = 0, end = 1, T0 = 1, x = 0.5, xi = 0.4;
            ASSERT(is_equal(0.2 , G(x, xi, T0, begin, end)));

            begin = 0, end = M_PI, T0 = 1, x = 0.4, xi = 0.5;
            ASSERT(std::abs(0.336338022763 - G(x, xi, T0, begin, end)) < eps);

            begin = 0, end = M_PI, T0 = 1, x = 0.5, xi = 0.4;
            ASSERT(std::abs(0.336338022763 - G(x, xi, T0, begin, end)) < eps);
        }

        void TestH(){

            ldouble eps = 1e-10;
            {
                ldouble begin = 0, end = 1,
                       T0 = 1, x = 0.5, xi = 0.25;
                std::vector<ldouble> t = linspace(begin, end, 1000000);
                ASSERT((std::abs(H(xi, x, t ,T0, begin, end) - H(xi, x, T0, begin, end)) < eps));
            }

            {
                ldouble begin = 0, end = 1,
                       T0 = 1, x = 0.25, xi = 0.5;
                std::vector<ldouble> t = linspace(begin, end, 1000000);
                ASSERT((std::abs(H(xi, x, t ,T0, begin, end) - H(xi, x, T0, begin, end)) < eps));
            }

            {
                ldouble begin = 0, end = 1,
                       T0 = 1, x = 0.5, xi = 0.25;
                size_t size = 1000000;
                std::vector<ldouble> t = linspace(begin, end, size);
                ldouble h = (end - begin)/(size - 1);
                t.front() += h/2;
                t.back() -= h/2;
                ASSERT((std::abs(H(xi, x, t ,T0, begin, end) - H(xi, x, T0, begin, end)) < eps));
            }

            {
                auto H_test = [](ldouble t, ldouble x, ldouble T0, ldouble end){
                    auto& l= end;
                    ldouble divider = T0*T0*l;
                    if(t < x){
                        return t*(l-x)*(2.0*l*x - t*t-x*x)/(6.0*divider);
                    }else{
                        return x*(l-t)*(2.0*l*t - t*t-x*x)/(6.0*divider);
                    }
                };

                ldouble begin = 0, end = 1, T0 = 1;
                std::vector<ldouble> t = linspace(begin, end, 1000);
                for(size_t i = 0; i < t.size(); ++i){
                    for(size_t j = 0; j < t.size(); ++j){
                        ASSERT((std::abs(H(t[j], t[i],T0, begin, end) - H_test(t[j], t[i], T0, end)) < eps));
                    }
                }
            }
        }

    }

    namespace membrane_dstrb {
        using namespace ::membrane_dstrb;

        bool CheckGrinFunction(ldouble x, ldouble y, ldouble xi, ldouble eta,
                               ldouble end_x, ldouble end_y, size_t M, size_t N,
                               ldouble eps, ldouble eps_loss, ldouble test_value){

            bool check = true;
            auto res = G(xi, eta, x, y, end_x, end_y, M, N);
            if(std::abs(res - test_value) >= eps_loss){
                std::cout << "G(xi, eta, x, y, end_x, end_y, M, N) = " << res << '\n'
                          << " loss = " << std::abs(res - test_value) << '\n';
                check = false;
            }
            res = G(xi, eta, x, y, end_x, end_y, eps);
            if(std::abs(res - test_value) >= eps_loss){
                std::cout << "G(xi, eta, x, y, end_x, end_y, eps) = " << res << '\n'
                          << " loss = " << std::abs(res - test_value) << '\n';
                check = false;
            }
            res = G(xi, eta, x, y, end_x, end_y, eps) - G(xi, eta, x, y, end_x, end_y, M, N);
            if(std::abs(res) >= eps_loss){
                std::cout << "G_eps - G_MN = " << res << '\n';
                check = false;
            }

            return check;
        }

        bool CheckHFunction(const Coordinates<ldouble> p, ldouble x, ldouble y, ldouble xi, ldouble eta,  ldouble begin_x,
                               ldouble end_x, ldouble begin_y, ldouble end_y, size_t M, size_t N,
                               ldouble eps, ldouble eps_loss, ldouble test_value){

            bool check = true;
            auto res = H(xi, eta, x, y, M, N, begin_x, end_x, begin_y, end_y);
            if(std::abs(res - test_value) >= eps_loss){
                std::cout << "H(xi, eta, x, y, M, N, begin_x, end_x, begin_y, end_y) = " << res << '\n'
                          << " loss = " << std::abs(res - test_value) << '\n';
                check = false;
            }
            res = H(xi, eta, x, y, eps, begin_x, end_x, begin_y, end_y);
            if(std::abs(res - test_value) >= eps_loss){
                std::cout << "H(xi, eta, x, y, eps, begin_x, end_x, begin_y, end_y) - test_value) = " << res << '\n'
                          << " loss = " << std::abs(res - test_value) << '\n';
                check = false;
            }
            res = H(xi, eta, x, y, M, N, begin_x, end_x, begin_y, end_y) -
                    H(p, xi, eta, x, y, M, N, begin_x, end_x, begin_y, end_y);
            if(std::abs(res) >= eps_loss){
                std::cout << "H_explicit - H_integral = " << res << '\n';
                check = false;
            }

            return check;
        }

        void TestGrinFunction(){
            {
                ldouble end_x = 1, end_y = 1;
                ldouble xi = 0.1, eta = 0.1, x = 0.1, y = 0.1,
                        eps = 1e-20, eps_loss = 1e-13, test_value = 0.0246709781804;
                size_t M = 10, N = 10;
                ASSERT(CheckGrinFunction(x, y, xi, eta, end_x, end_y, M, N, eps, eps_loss, test_value));
            }

            {
                 ldouble end_x = 1, end_y = 1;
                 ldouble xi = 0.4, eta = 0.3, x = 0.2, y = 0.7,
                         eps = 1e-30, eps_loss = 1e-5, test_value = 0.0778111431327;
                 size_t M = 10, N = 10;
                 ASSERT(CheckGrinFunction(x, y, xi, eta, end_x, end_y, M, N, eps, eps_loss, test_value));
            }


            {
                 ldouble end_x = 1, end_y = 2;
                 ldouble xi = 0.4, eta = 0.3, x = 0.2, y = 0.7,
                         eps = 1e-30, eps_loss = 1e-4, test_value = 0.215906582652;
                 size_t M = 10, N = 10;
                 ASSERT(CheckGrinFunction(x, y, xi, eta, end_x, end_y, M, N, eps, eps_loss, test_value));
            }


            {
                 ldouble end_x = 2, end_y = 1;
                 ldouble xi = 0.4, eta = 0.3, x = 0.2, y = 0.7,
                         eps = 1e-30, eps_loss = 1e-5, test_value = 0.174953519804;
                 size_t M = 10, N = 10;
                 ASSERT(CheckGrinFunction(x, y, xi, eta, end_x, end_y, M, N, eps, eps_loss, test_value));
            }
        }

        void TestH(){
            {
                size_t size_x = 100, size_y = size_x;
                ldouble begin_x = 0, end_x = 1.0,
                        begin_y = 0, end_y = 1.0, eps = 1e-20, eps_loss = 1e-14;
                ldouble xi = 0.1, eta = 0.1, x = 0.1, y = 0.1, test_value = 0.00018784499557;
                size_t M = 10, N = 10;
                auto x_vec = linspace(begin_x, end_x, size_x);
                auto y_vec = linspace(begin_y, end_y, size_y);
                ASSERT(CheckHFunction({x_vec, y_vec}, x, y, xi, eta,begin_x, end_x, begin_y, end_y, M, N, eps, eps_loss, test_value));
            }

            {
                size_t size_x = 100, size_y = size_x;
                ldouble begin_x = 0, end_x = 1.0,
                        begin_y = 0, end_y = 1.0, eps = 1e-30, eps_loss = 1e-5;
                ldouble xi = 0.4, eta = 0.3, x = 0.2, y = 0.7, test_value = 0.00562635387364;
                size_t M = 10, N = 10;
                auto x_vec = linspace(begin_x, end_x, size_x);
                auto y_vec = linspace(begin_y, end_y, size_y);
                ASSERT(CheckHFunction({x_vec, y_vec}, x, y, xi, eta,begin_x, end_x, begin_y, end_y, M, N, eps, eps_loss, test_value));
            }

            {
                size_t size_x = 100, size_y = size_x;
                ldouble begin_x = 0, end_x = 1.0,
                        begin_y = 0, end_y = 2.0, eps = 1e-30, eps_loss = 1e-7;
                ldouble xi = 0.4, eta = 0.3, x = 0.2, y = 0.7, test_value = 0.0574195957367;
                size_t M = 10, N = 10;
                auto x_vec = linspace(begin_x, end_x, size_x);
                auto y_vec = linspace(begin_y, end_y, size_y);
                ASSERT(CheckHFunction({x_vec, y_vec}, x, y, xi, eta,begin_x, end_x, begin_y, end_y, M, N, eps, eps_loss, test_value));
            }

            {
                size_t size_x = 100, size_y = size_x;
                ldouble begin_x = 0, end_x = 2.0,
                        begin_y = 0, end_y = 1.0, eps = 1e-30, eps_loss = 1e-12;
                ldouble xi = 0.4, eta = 0.3, x = 0.2, y = 0.7, test_value = 0.0373761237253;
                size_t M = 10, N = 10;
                auto x_vec = linspace(begin_x, end_x, size_x);
                auto y_vec = linspace(begin_y, end_y, size_y);
                ASSERT(CheckHFunction({x_vec, y_vec}, x, y, xi, eta,begin_x, end_x, begin_y, end_y, M, N, eps, eps_loss, test_value));
            }
        }

        void TestCalculateLoss(){
            {
                size_t size_x = 100, size_y = size_x;
                ldouble begin_x = 0, end_x = 1.0,
                        begin_y = 0, end_y = 1.0;

                auto x = linspace(begin_x, end_x, size_x);
                auto y = linspace(begin_y, end_y, size_y);

                std::vector<std::vector<ldouble>> u0(size_y, std::vector<ldouble>(size_x)),
                                                  u(size_y, std::vector<ldouble>(size_x));

                auto u0_function = [end_x, end_y](ldouble x, ldouble y){ return (x*(end_x-x))*y*(end_y-y)/10.0;};
                for(size_t i = 0 ; i < size_y; ++i){
                    for(size_t j = 0; j < size_x; ++j){
                        u0[i][j] = u0_function(x[j], y[i]);
                    }
                }
                u = u0;
                ASSERT_EQUAL(CalculateLoss({x, y}, u, u0), 0.0);
            }
        }
    }


    void TestGauseSolving(){
        {
            ldouble eps = 1e-16;
            std::vector<std::vector<ldouble>> matrix = { {299, 124, 167, 12789},
                                                        {3127, 125, 1274, 325},
                                                        {126, 1674, 162, 1678},
                                                        {869, 459, 1689, 214}
                                                      };
            std::vector<ldouble> y = { 125, 126, 175, 314 };
            SolveGauss(matrix, y);
            auto solution = y;
            std::vector<ldouble> solution_test = {-0.0379818975131009 , 0.0822681424354419 , 0.1821451807833707, 0.0074858935416060};

            for(size_t i = 0 ; i < solution.size(); ++i){
                ASSERT(std::abs(solution[i] - solution_test[i]) < eps)
            }
        }

        {
            size_t size = 100;
            std::vector<ldouble> b(size);
            std::ifstream fb("tests/b_matlab.csv");
            if(fb.is_open()){
                std::string line;
                size_t i = 0;
                while(std::getline(fb, line)){
                    b[i] = std::stold(line);
                    ++i;
                }
            }
            fb.close();

            std::vector<std::vector<ldouble>> matrix(size);
            std::string delimiter = ",";
            std::ifstream fmatrix("tests/matrix_matlab.csv");
            if(fmatrix.is_open()){
                std::string line;
                size_t i = 0;
                while(std::getline(fmatrix, line)){
                    matrix[i] = detail::SplitIntoNums(line, delimiter);
                    ++i;
                }
            }
            fmatrix.close();

            std::vector<ldouble> solution_matlab(size);
            std::ifstream fsolution("tests/solution_matlab.csv");
            if(fsolution.is_open()){
                std::string line;
                size_t i = 0;
                while(std::getline(fsolution, line)){
                    solution_matlab[i] = std::stold(line);
                    ++i;
                }
            }
            fsolution.close();

            ldouble eps = 1e-13;
            SolveGauss(matrix, b);
            auto solution = b;
            std::cout << std::setprecision(16);
            for(size_t i = 0; i < size; ++i){
                ldouble loss = std::abs(solution[i] - solution_matlab[i]);
                bool check = loss < eps;
                if(!check){
                    std::cout << "solution["<< i << "] = " << solution[i] << '\n';
                    std::cout << "solution_matlab["<< i << "] = " << solution_matlab[i] << '\n';
                    std::cout << "loss["<< i << "] = " << loss << '\n';
                }
                ASSERT(check);
            }

        }
    }

    void Test(){
        TestRunner runner;
        RUN_TEST(runner, TestIsEqual);
        RUN_TEST(runner, TestIsSegment);
        RUN_TEST(runner, TestLinspace);
        RUN_TEST(runner, str_dstrb::TestGrinFunction);
        RUN_TEST(runner, str_dstrb::TestH);
        RUN_TEST(runner, TestGauseSolving);
        RUN_TEST(runner, membrane_dstrb::TestCalculateLoss);
        RUN_TEST(runner, membrane_dstrb::TestGrinFunction);
        RUN_TEST(runner, membrane_dstrb::TestH);
    }
}


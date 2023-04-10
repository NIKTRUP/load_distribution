//#define NDEBUG

#include "include/test_framework.h"
#include "include/functions.h"
#include "include/tests.h"
#include "include/string_load_distribution.h"
#include <utility>
using namespace std;

// Записывает содержимое контейнера в csv файл
void Write(const membrane_dstrb::Coordinates<ldouble> p, std::vector<std::vector<ldouble>> u, const std::string& name){
    using namespace std::literals;
    std::ofstream out(name + ".csv"s);
    if(out.is_open()){
        for(size_t i = 0 ; i < p.y.size(); ++i){
            for(size_t j = 0 ; j < p.x.size(); ++j){
                out << p.x[j] << ',' << p.y[i]<< ',' << u[i][j] << '\n';
            }
        }
    }
    out.close();
}

void CheckLinspace(std::vector<ldouble>& x, ldouble begin, ldouble end, size_t size){
    ASSERT(end > begin);
    ASSERT(size > 2)
    ASSERT(x.size() == size)

    ldouble h = (end - begin)/(size-1);
    std::cout << std::setprecision(6);
    for(size_t i = 0 ; i < size - 1; ++i){
        ldouble eps = 1e-16;
        bool check_h = std::abs(x[i + 1] - x[i] - h) <= eps ;
        if(!check_h){
            std::cout << "abs(x["<< i << "] - " << " xi["<< i + 1 << "]) = " << std::abs(x[i] - x[i+1]) << '\n';
            std::cout << "h = " << h << '\n';
        }

        ASSERT(check_h);
    }
}

void StringDistributionExample(){
    // Init Data
    ldouble  begin = 0, end = 1, T0 = 1, P = 0.5, alpha = 100; // когда отношение P/T0 = 2, то отлично сходится
    size_t size = 51;
    auto u0_function = [end](ldouble x){ return (x*(end-x))/10.0; };
    auto x = linspace(begin, end, size);
    //ldouble h = (end - begin)/(size-1);
    // x.front() += h/10;
    // x.back() -= h/10;
    x = {x.begin() + 1, x.end() - 1};
    size = x.size();

    std::vector<std::vector<ldouble>> matrix(size + 1, std::vector<ldouble>(size+1));
    std::vector<ldouble> u0(size), u(size) ,b(matrix.size());

    // Fill Data
    for(size_t i = 0 ; i < size; ++i){
        u0[i] = u0_function(x[i]);
    }
    str_dstrb::FillSLAE(matrix, b, x, u0, alpha, P, T0, begin ,end);
    //WriteSLAE(matrix, b, "slae");

    // Solve Task
    SolveGauss(matrix, b);
    auto solution = b;

    auto p = std::vector<ldouble>{solution.begin(), solution.end() - 1};
    for(size_t i = 0; i < u.size(); ++i){
        u[i] = str_dstrb::СalculateDeviationFromEquilibrium(x[i], {x, p}, T0, begin, end);
    }
    std::cout << "loss = " << std::abs(str_dstrb::CalculateLoss(x, u , u0)) << std::endl;
    std::cout << "lambda = " << solution.back() << std::endl;

    x.push_back(end); x.insert(x.begin(), begin);
    u0.push_back(0); u0.insert(u0.begin(), 0);
    u.push_back(0); u.insert(u.begin(), 0);
    p.push_back(0); p.insert(p.begin(), 0);
    WriteLog(x, u0, "x_u0"s);
    WriteLog(x, u, "x_u"s);
    WriteLog(x, p, "x_p");
    std::cout << __FUNCTION__ << " is Done! \n" << std::endl;
}

void MembraneDistributionExample(){
    ldouble begin_x = 0, end_x = 1.0,
            begin_y = 0, end_y = 1.0,
            P = 0.5, alpha = 0.0;
    size_t size_x = 15, size_y = size_x,
            M = 5, N = 5;

    auto x = linspace(begin_x, end_x, size_x);
    x = {x.begin() + 1, x.end() - 1};
    size_x = x.size();

    auto y = linspace(begin_y, end_y, size_y);
    y = {y.begin() + 1, y.end() - 1};
    size_y = y.size();

    std::vector<std::vector<ldouble>> matrix(size_x*size_y+ 1, std::vector<ldouble>(size_x*size_y + 1, 0));
    std::vector<std::vector<ldouble>> u0(size_y, std::vector<ldouble>(size_x));
    std::vector<ldouble> b(matrix.size(), 0.0);

    auto u0_function = [end_x, end_y](ldouble x, ldouble y){ return (x*(end_x-x))*y*(end_y-y)/10.0;};
    for(size_t i = 0 ; i < size_y; ++i){
        for(size_t j = 0; j < size_x; ++j){
            u0[i][j] = u0_function(x[j], y[i]);
        }
    }

    membrane_dstrb::FillSLAE(matrix, b, {x, y}, u0, alpha, P, size_x, size_y, M, N, begin_x, end_x, begin_y, end_y);
    std::cout << "matrix is filled" << std::endl;
    WriteSLAE(matrix, b, "slae");

    SolveGauss(matrix, b);
    std::cout << "lambda = " << b.back() << std::endl;
    auto p = membrane_dstrb::FillP(b, size_x, size_y);
    std::cout << "slae is solved" << std::endl;

    auto u = membrane_dstrb::CalculateU({x, y, p}, M, N, end_x, end_y);
    x.push_back(end_x); x.insert(x.begin(), begin_x);
    y.push_back(end_y); y.insert(y.begin(), begin_y);
    u0 = membrane_dstrb::AddZerosBorders(u0);

    //std::cout << "loss = " << std::abs(membrane_dstrb::CalculateLoss({x, y}, u , u0)) << std::endl;

    p = membrane_dstrb::AddZerosBorders(p);
    Write({x, y}, u0, "u0");
    std::cout << "u0 is writed" << std::endl;
    Write({x, y}, u, "u");
    std::cout << "u is writed" << std::endl;
    Write({x, y}, p, "p");
    std::cout << "p is writed" << std::endl;
    std::cout << __FUNCTION__ << " is Done! \n" << std::endl;

}

int main(){

#ifndef NDEBUG
    tests::Test();
#endif
    //StringDistributionExample();
    //MembraneDistributionExample();
    std::cout << __FUNCTION__ << " is finished! \n" << std::endl;
}

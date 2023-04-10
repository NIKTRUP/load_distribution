#include "include/string_load_distribution.h"

namespace str_dstrb {

    // T0 - натяжение ненагруженной струны
    // Функция Грина
    ldouble G(ldouble x, ldouble xi, ldouble T0 ,ldouble begin, ldouble end){

        ldouble& l = end;
        if(is_equal(T0, 0.0) || is_equal(l, begin)){
            throw ArgumentError(" Невалидные аргументы"s);
        }

        if(in_segment(x, begin, xi)){
            return x*(l - xi)/(T0*l);
        }else if(in_segment(x, xi, l)){
            return xi*(l - x)/(T0*l);
        }else{
            throw ArgumentError(" x и xi находятся не в отрезке [begin, end]"s);
        }
    }

    // Функция H
    ldouble H(ldouble xi, ldouble x, std::vector<ldouble> t ,ldouble T0 ,ldouble begin, ldouble end){
        ldouble sum = 0;

        if(!in_segment(xi, begin, end)){
            throw ArgumentError("xi = " + std::to_string(xi) + " находится не в интервале от " + std::to_string(begin) + " до " + std::to_string(end));
        }else if (!in_segment(x, begin, end)){
            throw ArgumentError("x = " + std::to_string(x) + " находится не в интервале от " + std::to_string(begin) + " до " + std::to_string(end));
        }

        ldouble h = 0.0;
        for(size_t i = 0 ; i < t.size(); ++i){
            h = GetStep(t, i);

            sum += G(x, t[i], T0, begin,end)*G(t[i], xi, T0, begin, end)*h/2.0;

        }
        return sum;
    }

    // Функция H Явная
    ldouble H(ldouble xi, ldouble x, ldouble T0 , ldouble begin, ldouble end){

        if(!in_segment(xi, begin, end)){
            throw ArgumentError("xi = " + std::to_string(xi) + " находится не в интервале от " + std::to_string(begin) + " до " + std::to_string(end));
        }else if (!in_segment(x, begin, end)){
            throw ArgumentError("x = " + std::to_string(x) + " находится не в интервале от " + std::to_string(begin) + " до " + std::to_string(end));
        }

        ldouble& l = end;
        ldouble divider = T0*T0*l*l;
        if(xi <= x){
            return ((l-xi)*(l-x)*xi*xi*xi/3.0 +
                    (l-x)*xi*(l*(x*x - xi*xi)/2.0 - (x*x*x - xi*xi*xi)/3.0) +
                    xi*x*((l*l*l - x*x*x)/3.0 - l*l*x + l*x*x))/divider;
        }else if(x <= xi){
            return ((l-xi)*(l-x)*x*x*x/3.0 +
                    (l-xi)*x*(l*(xi*xi-x*x)/2.0 - (xi*xi*xi - x*x*x)/3.0) +
                    xi*x*((l*l*l - xi*xi*xi)/3.0 - l*l*xi + l*xi*xi))/divider;
        }else{
            throw CalculateError("Ошибка вычислений ");
        }
    }

    // TODO Tests
    // Вычисляет отклонение от положения равновесия
    ldouble СalculateDeviationFromEquilibrium(ldouble x, Points<ldouble> p,
                                             ldouble T0 , ldouble begin, ldouble end){
        size_t size = p.x.size();
        if(end <= begin || size <= 2|| p.x.size() != p.y.size()){
            throw ArgumentError(" Невалидные аргументы"s);
        }

        ldouble sum = 0.0;
        ldouble h = 0.0;
        h = p.x[1] - p.x[0];
        for(size_t i = 0; i < size; ++i){
            //h = GetStep(p.x, i);
            sum += G(x, p.x[i], T0, begin, end)*p.y[i]*h; // h/2.0
        }
        return sum;
    }

    // Вычисляет сумму квадрата ошибки между заданным отклонением и вычисленным
    ldouble CalculateLoss(std::vector<ldouble> x, std::vector<ldouble> u, std::vector<ldouble> u0){
        ldouble sum = 0.0, loss = 0.0, h = x[1] - x[0]; // h = 0.0;
        for(size_t i = 0 ; i < x.size(); ++i){
            //h = GetStep(x, i);

            loss = (u[i] - u0[i]);
            sum += (loss*loss)*h; // h/2.0
        }
        return sum;
    }

    // вычисляет расстояние между точками контейнера x
    ldouble GetStep(const std::vector<ldouble>& x, size_t i){
        ldouble h = 0.0;
        if(i == x.size() - 1){
            h = x[x.size() - 1] - x[x.size() - 2];
        }else if(i == 0){
            h = x[1] - x[0];
        }
        else{
            h = (x[i+1] - x[i - 1]);
        }
        return h;
    }

}

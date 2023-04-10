#include "include/functions.h"

#define EPS 1e-15

// проверяет на равенство два числа ldouble с точностью до 15 знаков после запятой
bool is_equal(ldouble x, ldouble y) {
    return std::fabs(x - y) < EPS;
}

// Проверяет принадлежность точки x к отрезку
bool in_segment(ldouble x, ldouble begin, ldouble end){
    return x - begin >= 0 && x - end <= EPS;
}

// Перегрузка
bool in_segment(ldouble x, std::pair<ldouble, ldouble> seg){
    return in_segment(x, seg.first, seg.second);
}

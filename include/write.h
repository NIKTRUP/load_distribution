#ifndef WRITE_H
#define WRITE_H
#include <string>
#include <fstream>
#include <iostream>
#include "include/membrane_load_distribution.h"

// сохраняет логи в csv формате
template <typename Container, typename Linspace>
void WriteLog(const Linspace& x, const Container& y, const std::string& name){
    using namespace std::literals;
    std::ofstream out(name + ".csv"s);
    if(out.is_open()){
        for(size_t i = 0 ; i < y.size(); ++i){
            out << x[i] << ',' << y[i] << '\n';
        }
    }
    out.close();
}

// Выводит логи в поток 'out'
template <typename Container, typename Linspace>
void WriteLog(const Linspace& x, const Container& y, std::ostream& out = std::cout){
    for(size_t i = 0 ; i < y.size(); ++i){
        out << x[i] << ',' << y[i] << '\n';
    }
}

// Записывает содержимое контейнера в csv файл
template <typename Container>
void Write(const Container& cont, const std::string& name){
    using namespace std::literals;
    std::ofstream out(name + ".csv"s);
    if(out.is_open()){
        out << cont << '\n';
    }
    out.close();
}

template <typename Matrix, typename Vector>
void WriteSLAE(const Matrix& matrix, const Vector& b, const std::string& name){
    using namespace std::literals;

    bool first = true;
    if(matrix.size() != b.size()){
        throw std::invalid_argument(" Количество строк и элементов правой части не совпадает ");
    }

    std::ofstream out(name + ".csv"s);
    if(out.is_open()){
        for(size_t i = 0; i < matrix.size(); ++i){
            if (first) {
                out << matrix[i] << ",=," << b[i];
                first = false;
                continue;
            }
            out << '\n' << matrix[i] << ",=," << b[i];
        }
    }
    out.close();
}

#endif // WRITE_H

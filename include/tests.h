#ifndef TESTS_H
#define TESTS_H

#include "include/test_framework.h"
#include "include/functions.h"
#include "include/string_load_distribution.h"
#include <iomanip>

namespace tests {
    namespace detail {

        void GenerateLinspace(ldouble begin, ldouble end, size_t size_dots);

        std::vector<ldouble> SplitIntoNums(std::string_view str, std::string delimiter = ",");

        std::tuple<std::vector<std::vector<ldouble>>, std::vector<ldouble>> GenerateSLAE(size_t size);
    }

    void TestLinspace();

    void TestIsEqual();

    void TestIsSegment();

    namespace str_dstrb {

        void TestGrinFunction();

        void TestH();

        void TestGetStep();
    }

    namespace membrane_dstrb {

        void TestGrinFunction();

        void TestH();

        void TestCalculateLoss();
    }

    void TestGauseSolving();

    void Test();
}

#endif // TESTS_H

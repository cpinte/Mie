#pragma once

#include <vector>
#include <complex>

namespace Mie {

    struct Result {
        std::vector<std::complex<double>> s1;
        std::vector<std::complex<double>> s2;
        double q_ext;
        double q_sca;
        double q_back;
        double g_sca;
    };

    Result bhmie(double x, const std::complex<double>& refrel, int nang);

}

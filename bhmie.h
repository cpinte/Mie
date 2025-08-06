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

    struct MuellerResult {
      double q_ext;
      double q_sca;
      double g_sca;
      std::vector<double> s11;
      std::vector<double> s12;
      std::vector<double> s22;
      std::vector<double> s33;
      std::vector<double> s34;
      std::vector<double> s44;
    };

    MuellerResult mueller_mie(double x, std::complex<double> refrel, const bool only_g);

}

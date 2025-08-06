// 1. Include guard to prevent multiple inclusions
#ifndef BHMIE_H
#define BHMIE_H

#include <vector>
#include <complex>

// 2. Namespace to prevent naming collisions
namespace Mie {

    // The result structure is now inside the Mie namespace
    struct Result {
        std::vector<std::complex<double>> s1;
        std::vector<std::complex<double>> s2;
        double q_ext;
        double q_sca;
        double q_back;
        double g_sca;
    };

    // The function declaration is also in the namespace
    Result bhmie(double x, const std::complex<double>& refrel, int nang);

} // namespace Mie

#endif // BHMIE_H

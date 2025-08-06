#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <numbers> // For std::numbers::pi

/**
 * @brief A structure to hold all the output results from the Mie calculation.
 */
struct MieResult {
    std::vector<std::complex<double>> s1;
    std::vector<std::complex<double>> s2;
    double q_ext;
    double q_sca;
    double q_back;
    double g_sca;
};

/**
 * @brief C++ translation of the Bohren-Huffman BHMIE Fortran subroutine.
 *
 * This function calculates scattering and absorption properties of a homogenous,
 * isotropic sphere using Mie theory.
 *
 * @param x The size parameter (2 * pi * radius / wavelength).
 * @param refrel The complex refractive index of the sphere relative to the surrounding medium.
 * @param nang The number of angles between 0 and 90 degrees for which to calculate scattering amplitudes.
 * @return A MieResult struct containing the scattering efficiencies and complex scattering amplitudes.
 */
MieResult bhmie(double x, const std::complex<double>& refrel, int nang) {
    // Ensure nang is at least 2 for calculations at 0, 90, and 180 degrees.
    if (nang < 2) {
        nang = 2;
    }

    // Local variables
    const std::complex<double> y = x * refrel;
    const double y_mod = std::abs(y);

    // Determine the number of terms to sum in the series
    const double x_stop_d = x + 4.0 * std::pow(x, 1.0/3.0) + 2.0;
    const int n_stop = static_cast<int>(x_stop_d);
    const int n_mx = static_cast<int>(std::max(x_stop_d, y_mod) + 15.0);

    // Dynamically sized vectors
    std::vector<std::complex<double>> d(n_mx + 1); // Fortran D(NMX) -> C++ d[n_mx]
    std::vector<double> amu(nang);
    std::vector<double> pi(nang);
    std::vector<double> pi0(nang, 0.0);
    std::vector<double> pi1(nang, 1.0);
    std::vector<double> tau(nang);

    // Initialize the result structure
    MieResult result;
    const int s_size = 2 * nang - 1;
    result.s1.assign(s_size, {0.0, 0.0});
    result.s2.assign(s_size, {0.0, 0.0});

    // Set up angles
    for (int j = 0; j < nang; ++j) {
        amu[j] = std::cos(static_cast<double>(j) * (std::numbers::pi / 2.0) / (nang - 1.0));
    }

    // --- Main Calculation ---

    // Logarithmic derivative D(n) calculated by downward recurrence
    d[n_mx] = {0.0, 0.0};
    for (int n = n_mx - 1; n >= 0; --n) {
        const double en = static_cast<double>(n + 1);
        d[n] = (en / y) - (1.0 / (d[n + 1] + en / y));
    }

    // Riccati-Bessel functions with real argument X calculated by upward recurrence
    double psi0 = std::cos(x);
    double psi1 = std::sin(x);
    double chi0 = -std::sin(x);
    double chi1 = std::cos(x);
    std::complex<double> xi1 = {psi1, -chi1};

    result.q_sca = 0.0;
    result.g_sca = 0.0;

    std::complex<double> an, bn, an1, bn1;
    double p = -1.0;

    for (int n = 1; n <= n_stop; ++n) {
        const double en = static_cast<double>(n);
        const double fn = (2.0 * en + 1.0) / (en * (en + 1.0));

        // Calculate psi_n and chi_n
        const double psi = (2.0 * en - 1.0) * psi1 / x - psi0;
        const double chi = (2.0 * en - 1.0) * chi1 / x - chi0;
        const std::complex<double> xi = {psi, -chi};

        // Store previous values of AN and BN for g_sca calculation
        if (n > 1) {
            an1 = an;
            bn1 = bn;
        }

        // Compute Mie coefficients AN and BN
        an = (d[n] / refrel + en / x) * psi - psi1;
        an /= ((d[n] / refrel + en / x) * xi - xi1);
        bn = (refrel * d[n] + en / x) * psi - psi1;
        bn /= ((refrel * d[n] + en / x) * xi - xi1);

        // Augment sums for Q_sca and g_sca
        result.q_sca += (2.0 * en + 1.0) * (std::norm(an) + std::norm(bn));
        result.g_sca += ((2.0 * en + 1.0) / (en * (en + 1.0))) * (an.real() * bn.real() + an.imag() * bn.imag());
        if (n > 1) {
            result.g_sca += ((en - 1.0) * (en + 1.0) / en) *
                           (an1.real() * an.real() + an1.imag() * an.imag() +
                            bn1.real() * bn.real() + bn1.imag() * bn.imag());
        }

        // Calculate scattering intensity pattern for angles 0 to 90
        for (int j = 0; j < nang; ++j) {
            pi[j] = pi1[j];
            tau[j] = en * amu[j] * pi[j] - (en + 1.0) * pi0[j];
            result.s1[j] += fn * (an * pi[j] + bn * tau[j]);
            result.s2[j] += fn * (an * tau[j] + bn * pi[j]);
        }

        // Use symmetry for angles > 90
        p = -p;
        for (int j = 0; j < nang - 1; ++j) {
            int jj = 2 * nang - 2 - j;
            result.s1[jj] += fn * p * (an * pi[j] - bn * tau[j]);
            result.s2[jj] += fn * p * (bn * pi[j] - an * tau[j]);
        }

        // Update recurrence relations
        psi0 = psi1;
        psi1 = psi;
        chi0 = chi1;
        chi1 = chi;
        xi1 = {psi1, -chi1};

        // Update pi_n for next value of n
        for (int j = 0; j < nang; ++j) {
            double temp_pi1 = ((2.0 * en + 1.0) * amu[j] * pi[j] - (en + 1.0) * pi0[j]) / en;
            pi0[j] = pi[j];
            pi1[j] = temp_pi1;
        }
    }

    // Finalize calculations
    result.g_sca = 2.0 * result.g_sca / result.q_sca;
    result.q_sca = (2.0 / (x * x)) * result.q_sca;
    result.q_ext = (4.0 / (x * x)) * result.s1[0].real();
    result.q_back = (4.0 / (x * x)) * std::norm(result.s1[2 * nang - 2]);

    return result;
}

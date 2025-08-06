#pragma once

#include <vector>
#include <complex>
#include <optional>

namespace Mie {

    /**
     * @brief Result of Mie scattering calculation.
     *
     * Contains scattering amplitudes and efficiencies.
     */
    struct Result {
        std::vector<std::complex<double>> s1;
        std::vector<std::complex<double>> s2;
        double q_ext;
        double q_sca;
        double q_back;
        double g_sca;
    };

    /**
     * @brief Computes Mie scattering for a sphere.
     *
     * @param x Size parameter (2 * pi * radius / wavelength).
     * @param refrel Complex refractive index of the sphere relative to the medium.
     * @param nang Number of angles for scattering amplitude calculation.
     * @return Result struct with scattering data.
     */
    Result bhmie(double x, const std::complex<double>& refrel, int nang);

    /**
     * @brief Mueller matrix result from Mie scattering.
     *
     * Contains Mueller matrix elements and efficiencies.
     */
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

    /**
     * @brief Computes the Mueller matrix for Mie scattering.
     *
     * @param x Size parameter (2 * pi * radius / wavelength).
     * @param refrel Complex refractive index of the sphere relative to the medium.
     * @param only_g If true, only compute asymmetry parameter g.
     * @return MuellerResult struct with matrix elements and efficiencies.
     */
    MuellerResult mueller_mie(double x, std::complex<double> refrel,  const std::optional<bool> only_g = false);

}

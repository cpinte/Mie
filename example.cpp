#include "bhmie.h"
#include <iostream>
#include <complex>
#include <numbers>


/**
 * @brief Example usage of the Mie scattering functions.
 *
 * Demonstrates how to use Mie::bhmie and Mie::mueller_mie to compute scattering properties
 * and print results for a sphere in a given medium.
 */
int main() {
    // --- Input Parameters ---
    const double radius_microns = 0.5;
    const double wavelength_microns = 0.6328; // HeNe laser
    const std::complex<double> refractive_index_sphere = {1.59, 0.0};
    const double refractive_index_medium = 1.0; // Vacuum or air
    const int num_angles = 91; // Calculate for 91 angles from 0 to 90 degrees

    // --- Derived Parameters ---
    const double x = 2.0 * std::numbers::pi * radius_microns / wavelength_microns;
    const std::complex<double> refrel = refractive_index_sphere / refractive_index_medium;

    try {

      // --- Call the Mie scattering function ---
      Mie::Result results;
      results  = Mie::bhmie(x, refrel, num_angles);

      // --- Print Results ---
      std::cout << "Mie Scattering Calculation Results" << std::endl;
      std::cout << "----------------------------------" << std::endl;
      std::cout << "Size Parameter (x): " << x << std::endl;
      std::cout << "Relative Refractive Index: " << refrel << std::endl;
      std::cout << std::endl;
      std::cout << "Efficiencies:" << std::endl;
      std::cout << "  Q_ext:  " << results.q_ext << std::endl;
      std::cout << "  Q_sca:  " << results.q_sca << std::endl;
      std::cout << "  Q_back: " << results.q_back << std::endl;
      std::cout << "  g_sca (Asymmetry Parameter): " << results.g_sca << std::endl;
      std::cout << std::endl;
      std::cout << "Scattering Amplitudes at key angles:" << std::endl;
      std::cout << "  Angle    S1                     S2" << std::endl;
      std::cout << "  -----    --                     --" << std::endl;
      std::cout << "  0 deg:   " << results.s1[0] << "    " << results.s2[0] << std::endl;
      std::cout << "  90 deg:  " << results.s1[num_angles - 1] << "    " << results.s2[num_angles - 1] << std::endl;
      std::cout << "  180 deg: " << results.s1[2 * num_angles - 2] << " " << results.s2[2 * num_angles - 2] << std::endl;


    } catch (const std::exception& e) {
      std::cerr << "BH Mie error: " << e.what() << std::endl;
    }


    try {
        // Call the new mueller_mie function
      Mie::MuellerResult mueller = Mie::mueller_mie(x, refrel);

        std::cout << "Mueller Matrix Calculation Results" << std::endl;
        std::cout << "----------------------------------" << std::endl;
        std::cout << "Q_ext: " << mueller.q_ext << std::endl;
        std::cout << "Q_sca: " << mueller.q_sca << std::endl;
        std::cout << "g_sca: " << mueller.g_sca << std::endl;
        std::cout << std::endl;

        if (!mueller.s11.empty()) {
            std::cout << "Mueller Matrix Elements at 0, 90, and 180 degrees:" << std::endl;
            std::cout << "Angle |    S11   |    S12   " << std::endl;
            std::cout << "----------------------------------" << std::endl;
            // Angle 0 degrees (index 0)
            std::cout << "  0   | " << mueller.s11[0] << " | " << mueller.s12[0] << std::endl;
            // Angle 90 degrees (index 90)
            std::cout << " 90   | " << mueller.s11[90] << " | " << mueller.s12[90] << std::endl;
            // Angle 180 degrees (index 180)
            std::cout << " 180  | " << mueller.s11[180] << " | " << mueller.s12[180] << std::endl;
        }


    } catch (const std::exception& e) {
      std::cerr << "Mueller error: " << e.what() << std::endl;
    }

    return 0;
}

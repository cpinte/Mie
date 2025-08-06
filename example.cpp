// --- Example Usage ---
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

    // --- Call the Mie scattering function ---
    MieResult results = bhmie(x, refrel, num_angles);

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

    return 0;
}

# Mie Scattering Library

A modern C++ implementation of the Bohren & Huffman Mie scattering algorithm with Python bindings. This library provides efficient computation of light scattering by spherical particles, including both amplitude scattering functions and Mueller matrix elements.

## Features

- **Mie Scattering Calculation**: Compute scattering amplitudes (S1, S2) and efficiencies (Q_ext, Q_sca, Q_back, g_sca)
- **Mueller Matrix Support**: Calculate full Mueller matrix elements for polarized light scattering
- **C++20 Implementation**: Modern C++ with comprehensive error handling
- **Python Bindings**: Easy-to-use Python interface via pybind11
- **Cross-platform**: Works on Windows, macOS, and Linux
- **MIT Licensed**: Free for academic and commercial use

## Installation

### Prerequisites

- CMake 3.16 or higher
- C++20 compatible compiler
- Python 3.7+ abd pybind11


### Building from Source

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd Mie
   ```

2. **Create build directory and configure**:
   ```bash
   mkdir build
   cd build
   cmake ..
   ```

3. **Build the project**:
   ```bash
   make
   ```

This will create:
- `mie_example` - C++ executable demonstrating the library usage
- `bhmie_py.*.so` (or `.pyd` on Windows) - Python module

## Usage

### C++ Example

```cpp
#include "bhmie.h"
#include <complex>
#include <numbers>

int main() {
    // Parameters
    const double radius_microns = 0.5;
    const double wavelength_microns = 0.6328; // HeNe laser
    const std::complex<double> refractive_index_sphere = {1.59, 0.0};
    const double refractive_index_medium = 1.0;
    const int num_angles = 91;

    // Calculate size parameter
    const double x = 2.0 * std::numbers::pi * radius_microns / wavelength_microns;
    const std::complex<double> refrel = refractive_index_sphere / refractive_index_medium;

    // Compute Mie scattering
    Mie::Result results = Mie::bhmie(x, refrel, num_angles);

    // Access results
    std::cout << "Q_ext: " << results.q_ext << std::endl;
    std::cout << "Q_sca: " << results.q_sca << std::endl;
    std::cout << "g_sca: " << results.g_sca << std::endl;

    // Compute Mueller matrix
    Mie::MuellerResult mueller = Mie::mueller_mie(x, refrel);
    std::cout << "S11 at 0°: " << mueller.s11[0] << std::endl;

    return 0;
}
```

### Python Example

```python
import bhmie_py
import numpy as np

# Compute Mie scattering
res = bhmie_py.bhmie(2.0, 1.5 + 0.1j, 90)
print(f"Scattering efficiency: {res.q_sca}")

# Compute Mueller matrix
mueller = bhmie_py.mueller_mie(2.0, 1.5 + 0.1j)
print(f"S11 elements: {mueller.s11[:5]}")  # First 5 elements
```

## API Reference

### C++ API

#### `Mie::Result` Structure
```cpp
struct Result {
    std::vector<std::complex<double>> s1;  // S1 scattering amplitude
    std::vector<std::complex<double>> s2;  // S2 scattering amplitude
    double q_ext;                          // Extinction efficiency
    double q_sca;                          // Scattering efficiency
    double q_back;                         // Backscattering efficiency
    double g_sca;                          // Asymmetry parameter
};
```

#### `Mie::MuellerResult` Structure
```cpp
struct MuellerResult {
    double q_ext;                          // Extinction efficiency
    double q_sca;                          // Scattering efficiency
    double g_sca;                          // Asymmetry parameter
    std::vector<double> s11;               // Mueller matrix element S11
    std::vector<double> s12;               // Mueller matrix element S12
    std::vector<double> s22;               // Mueller matrix element S22
    std::vector<double> s33;               // Mueller matrix element S33
    std::vector<double> s34;               // Mueller matrix element S34
    std::vector<double> s44;               // Mueller matrix element S44
};
```

#### Functions

**`Mie::bhmie(x, refrel, nang)`**
- **Parameters**:
  - `x`: Size parameter (2π × radius / wavelength)
  - `refrel`: Complex refractive index of sphere relative to medium
  - `nang`: Number of angles for scattering amplitude calculation
- **Returns**: `Mie::Result` with scattering data

**`Mie::mueller_mie(x, refrel, only_g)`**
- **Parameters**:
  - `x`: Size parameter
  - `refrel`: Complex refractive index
  - `only_g`: Optional boolean to compute only asymmetry parameter
- **Returns**: `Mie::MuellerResult` with Mueller matrix elements

### Python API

The Python module `bhmie_py` provides the same functionality with identical function signatures:

```python
# Mie scattering
result = bhmie_py.bhmie(x, refrel, nang)

# Mueller matrix
mueller = bhmie_py.mueller_mie(x, refrel, only_g=None)
```


## Examples

Run the included examples:

```bash
# C++ example
./mie_example

# Python example
python example.py
```

## Documentation

Detailed API documentation is available in the `docs/` directory. To generate documentation:

```bash
doxygen Doxyfile
```

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

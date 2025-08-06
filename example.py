"""
Example usage of the bhmie_py Python bindings for Mie scattering.

This script demonstrates how to compute Mie scattering and the Mueller matrix
using the bhmie_py module, which provides Python bindings to the C++ implementation.
"""
import bhmie_py

# Compute Mie scattering for a sphere with size parameter 2.0 and complex refractive index 1.5 + 0.1j
res = bhmie_py.bhmie(2.0, 1.5 + 0.1j, 90)
print(res.q_sca)

# Compute the Mueller matrix for the same parameters
mueller = bhmie_py.mueller_mie(2.0, 1.5 + 0.1j)
print(mueller.s11)

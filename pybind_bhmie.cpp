#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <optional>
#include "bhmie.h"

namespace py = pybind11;

/**
 * @brief Python bindings for the BHMIE Mie scattering module using pybind11.
 *
 * Exposes the Result and MuellerResult structs and the bhmie and mueller_mie functions to Python.
 */
PYBIND11_MODULE(bhmie_py, m) {
    m.doc() = "Python bindings for the BHMIE Mie scattering module";

    // Bind Result struct
    py::class_<Mie::Result>(m, "Result")
        .def_readonly("s1", &Mie::Result::s1)
        .def_readonly("s2", &Mie::Result::s2)
        .def_readonly("q_ext", &Mie::Result::q_ext)
        .def_readonly("q_sca", &Mie::Result::q_sca)
        .def_readonly("q_back", &Mie::Result::q_back)
        .def_readonly("g_sca", &Mie::Result::g_sca);

    // Bind MuellerResult struct
    py::class_<Mie::MuellerResult>(m, "MuellerResult")
        .def_readonly("q_ext", &Mie::MuellerResult::q_ext)
        .def_readonly("q_sca", &Mie::MuellerResult::q_sca)
        .def_readonly("g_sca", &Mie::MuellerResult::g_sca)
        .def_readonly("s11", &Mie::MuellerResult::s11)
        .def_readonly("s12", &Mie::MuellerResult::s12)
        .def_readonly("s22", &Mie::MuellerResult::s22)
        .def_readonly("s33", &Mie::MuellerResult::s33)
        .def_readonly("s34", &Mie::MuellerResult::s34)
        .def_readonly("s44", &Mie::MuellerResult::s44);

    // Bind functions
    m.def("bhmie", &Mie::bhmie,
          py::arg("x"), py::arg("refrel"), py::arg("nang"),
          "Compute Mie scattering using the BHMIE algorithm");

    m.def("mueller_mie", &Mie::mueller_mie,
          py::arg("x"), py::arg("refrel"), py::arg("only_g") = std::nullopt,
          "Compute Mueller matrix from Mie scattering");
}

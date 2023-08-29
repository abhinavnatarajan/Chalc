#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#include <pybind11/pybind11.h>

namespace
{
    namespace py = pybind11;
}

PYBIND11_MODULE(_version, m) {
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
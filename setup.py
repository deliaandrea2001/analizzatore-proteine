# setup.py - compila il core C++ "seq_core".
#   python setup.py build_ext --inplace

import sys
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Su Linux il link di un'estensione C++ a volte avviene con gcc senza la
# runtime C++: forziamo -lstdc++ per evitare "undefined symbol _ZTVN10...".
# Su macOS NON serve (usa libc++) e romperebbe il build, quindi lo escludiamo.
extra_link = ["-lstdc++"] if sys.platform.startswith("linux") else []

ext_modules = [
    Pybind11Extension("seq_core", ["seq_core.cpp"], cxx_std=14,
                      extra_link_args=extra_link),
]

setup(
    name="seq_core",
    version="0.1",
    description="Core C++ per analisi di sequenze proteiche",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)

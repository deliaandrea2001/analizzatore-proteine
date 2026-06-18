# setup.py - compila il core C++ "seq_core".
#   python setup.py build_ext --inplace

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension("seq_core", ["seq_core.cpp"], cxx_std=14),
]

setup(
    name="seq_core",
    version="0.1",
    description="Core C++ per analisi di sequenze proteiche",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)

from setuptools import setup, Extension
import setuptools_scm  # noqa  Ensure it’s installed
from setuptools_rust import Binding, RustExtension

extensions = [
    Extension("cutadapt._align", sources=["src/cutadapt/_align.pyx"]),
]
rust_extensions = [
    RustExtension("cutadapt.qualtrim", "src/rust/Cargo.toml", binding=Binding.PyO3)
]
setup(ext_modules=extensions, rust_extensions=rust_extensions)

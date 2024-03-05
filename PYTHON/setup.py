import torch
import setuptools
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

setup(
    name='mfdfa',
    ext_modules=[
        CUDAExtension(
            name='mfdfa',
            sources=['mfdfa.cpp'],
            extra_compile_args={'cxx':['-O3'],})
    ],
    cmdclass={
        'build_ext': BuildExtension
})

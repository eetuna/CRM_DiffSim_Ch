"""Setup script for crm_ml_rl package."""

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import sys
import os
import subprocess


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake required")
        
        for ext in self.extensions:
            self.build_extension(ext)
    
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        build_dir = os.path.join('crm_ml_rl', 'bindings', 'build')
        os.makedirs(build_dir, exist_ok=True)
        
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            '-DCMAKE_BUILD_TYPE=Release'
        ]
        
        build_args = ['--', '-j8']
        
        subprocess.check_call(['cmake', os.path.join('..', '..')] + cmake_args, cwd=build_dir)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=build_dir)


setup(
    name='crm_ml_rl',
    version='1.0.0',
    author='Your Name',
    description='ML/RL for MRI-Actuated Robotic Catheter',
    packages=find_packages(),
    ext_modules=[CMakeExtension('crm_cpp', sourcedir='crm_ml_rl/bindings')],
    cmdclass={'build_ext': CMakeBuild},
    install_requires=[
        'numpy<2.0.0',
        'torch>=2.0.0',
        'gymnasium>=0.29.0',
        'stable-baselines3>=2.2.0',
        'pybind11>=2.11.0',
    ],
    python_requires='>=3.10',
    zip_safe=False,
)

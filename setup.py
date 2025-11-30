"""Setup configuration for CRM Catheter Python package."""

import os
import sys
import subprocess
from pathlib import Path
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """Extension that uses CMake to build."""
    
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """Build extension using CMake."""
    
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        
        cmake_args = [
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
            f'-DPYTHON_EXECUTABLE={sys.executable}',
            '-DBUILD_PYTHON_BINDINGS=ON',
            '-DBUILD_TESTING=OFF',
            '-DBUILD_EXAMPLES=OFF',
        ]
        
        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        cmake_args += [f'-DCMAKE_BUILD_TYPE={cfg}']
        
        if 'CMAKE_BUILD_PARALLEL_LEVEL' not in os.environ:
            build_args += ['-j4']
        
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        
        subprocess.check_call(
            ['cmake', ext.sourcedir] + cmake_args,
            cwd=self.build_temp
        )
        
        subprocess.check_call(
            ['cmake', '--build', '.', '--target', 'crm_cpp'] + build_args,
            cwd=self.build_temp
        )


long_description = (Path(__file__).parent / "README.md").read_text(encoding="utf-8")
requirements = (Path(__file__).parent / "requirements.txt").read_text().splitlines()
requirements = [r.strip() for r in requirements if r.strip() and not r.startswith('#')]

setup(
    name='crm-catheter',
    version='1.0.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='Cosserat Rod Model for MRI-Actuated Robotic Catheter with ML/AI Integration',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/CRM_Dynamics',
    project_urls={
        'Bug Tracker': 'https://github.com/yourusername/CRM_Dynamics/issues',
        'Documentation': 'https://crm-catheter.readthedocs.io/',
        'Source Code': 'https://github.com/yourusername/CRM_Dynamics',
    },
    packages=find_packages(where='python', exclude=['tests', 'examples']),
    package_dir={'': 'python'},
    ext_modules=[CMakeExtension('crm_cpp')],
    cmdclass={'build_ext': CMakeBuild},
    install_requires=requirements,
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    keywords='robotics continuum catheter machine-learning reinforcement-learning',
    zip_safe=False,
)

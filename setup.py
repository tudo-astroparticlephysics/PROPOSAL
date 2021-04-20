# see https://stackoverflow.com/a/48015772/3838691
import os
import sys
import sysconfig
import pathlib
import subprocess as sp
import re

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


SETUP_DIR = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(SETUP_DIR, 'CMakeLists.txt')) as f:
    f_read = f.read()
    res = re.search(r'PROPOSAL_VERSION_MAJOR (\d)', f_read)
    v_major = res.group(1)
    res = re.search(r'PROPOSAL_VERSION_MINOR (\d)', f_read)
    v_minor = res.group(1)
    res = re.search(r'PROPOSAL_VERSION_PATCH (\d)', f_read)
    v_patch = res.group(1)
    version = v_major + '.' + v_minor + '.' + v_patch


def get_cmake():
    ''' On centos 7, cmake is cmake 2.x but we need > 3.8'''
    for exe in ['cmake3', 'cmake']:
        try:
            ret = sp.run([exe, '--version'], stdout=sp.PIPE, stderr=sp.PIPE)
        except Exception:
            continue

        if ret.returncode == 0:
            return exe
    raise OSError('You need cmake >= 3.9')


class CMakeExtension(Extension):
    def __init__(self, name, source_dir=None, target=None, **kwargs):
        if source_dir is None:
            self.source_dir = SETUP_DIR
        else:
            self.source_dir = os.path.join(SETUP_DIR, source_dir)
        self.target = target
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[], **kwargs)


class build_ext_cmake(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)

    def build_cmake(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        for d in (self.build_temp, extdir):
            os.makedirs(d, exist_ok=True)

        cfg = 'Debug' if self.debug else 'Release'
        cmake = get_cmake()

        rpath = '@loader_path' if sys.platform == 'darwin' else '$ORIGIN'
        python_lib = pathlib.Path(
            sysconfig.get_config_var('LIBDIR'),
            sysconfig.get_config_var('INSTSONAME')
        )
        if not os.getenv('NO_CONAN', False):
            print("Using conan to install dependencies. Set environment variable NO_CONAN to skip conan.")
            conan_call = [
                'conan',
                'install',
                ext.source_dir,
                '-o with_python=True',
                '-o with_testing=False'
            ]
            sp.run(conan_call, cwd=self.build_temp, check=True)
        cmake_call = [
            cmake,
            ext.source_dir,
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            '-DCMAKE_BUILD_TYPE=' + cfg,
            '-DBUILD_TESTING=OFF',
            '-DBUILD_PYTHON=ON',
            '-DCMAKE_POSITION_INDEPENDENT_CODE=TRUE',
            '-DBUILD_EXAMPLE=OFF',
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DPYTHON_LIBRARY=' + str(python_lib),
            '-DCMAKE_INSTALL_RPATH={}'.format(rpath),
            '-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON',
            '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=OFF',
            '-DPYTHON_INCLUDE_DIR=' + sysconfig.get_path('include'),
        ]
        sp.run(cmake_call, cwd=self.build_temp, check=True)
        build_call = [
            cmake,
            '--build', '.',
            '--config', cfg,
            '--', '-j{}'.format(os.getenv('BUILD_CORES', 2))
        ]
        if ext.target is not None:
            build_call.extend(['--target', ext.target])
        sp.run(build_call, cwd=self.build_temp, check=True)


setup(
    version=version,
    ext_modules=[
        CMakeExtension('proposal'),
    ],
    cmdclass={
        'build_ext': build_ext_cmake,
    },
    zip_safe=False,
)

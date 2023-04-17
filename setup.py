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
with open(os.path.join(SETUP_DIR, "CMakeLists.txt")) as f:
    f_read = f.read()
    res = re.search(r"PROPOSAL_VERSION_MAJOR (\d)", f_read)
    v_major = res.group(1)
    res = re.search(r"PROPOSAL_VERSION_MINOR (\d)", f_read)
    v_minor = res.group(1)
    res = re.search(r"PROPOSAL_VERSION_PATCH (\d)", f_read)
    v_patch = res.group(1)
    version = v_major + "." + v_minor + "." + v_patch


def get_cmake():
    """ On centos 7, cmake is cmake 2.x but we need > 3.8"""
    for exe in ["cmake3", "cmake"]:
        try:
            ret = sp.run([exe, "--version"], stdout=sp.PIPE, stderr=sp.PIPE)
        except Exception:
            continue

        if ret.returncode == 0:
            return exe
    raise OSError("You need cmake >= 3.16")

def exists_conan_default_file():
    profiles = sp.check_output(["conan", "profile", "list"], encoding="UTF-8").split()
    if "default" in profiles:
        return True
    return False

def create_conan_profile():
    cmd = ["conan", "profile", "detect"]
    r = sp.run(cmd)
    if r.returncode != 0:
        raise RuntimeError(
            "conan was not able to create a defaultp profile."
        )

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
        CMAKE_CXX_FLAGS = ""
        if not os.getenv("NO_CONAN", False):
            print(
                "Using conan to install dependencies. Set environment variable NO_CONAN to skip conan."
            )

            if not exists_conan_default_file():
                create_conan_profile()

            conan_call = [
                'conan',
                'build',
                ext.source_dir,
                '-o with_python=True',
                '-o with_testing=False',
                '--build=missing'
               ]
            sp.run(conan_call, cwd=self.build_temp, check=True)

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

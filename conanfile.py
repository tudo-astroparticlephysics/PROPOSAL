from conans import ConanFile, CMake, tools
import os


class PROPOSALConan(ConanFile):
    name = "proposal"
    homepage = "https://github.com/tudo-astroparticlephysics/PROPOSAL"
    license = "GNU Lesser General Public License v3.0"
    description = "the very best lepton and photon propagator"
    topics = ("propagator", "lepton", "photon", "stochastic")
    settings = "os", "compiler", "build_type", "arch"
    exports_sources = "*"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "with_testing": [True, False],
        "with_python": [True, False],
        "with_documentation": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_testing": False,
        "with_python": False,
        "with_documentation": False,
    }
    generators = "cmake_find_package", "cmake_paths"
    _cmake = None

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def requirements(self):
        self.requires("cubicinterpolation/0.1.5")
        self.requires("spdlog/1.10.0")
        self.requires("nlohmann_json/3.9.1")
        if self.options.with_python:
            self.requires("pybind11/2.6.2")
        if self.options.with_testing:
            self.requires("boost/1.78.0")
            self.requires("gtest/1.11.0")
        if self.options.with_documentation:
            self.requires("doxygen/1.8.20")

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake
        self._cmake = CMake(self)
        self._cmake.definitions["BUILD_TESTING"] = self.options.with_testing
        self._cmake.definitions["BUILD_PYTHON"] = self.options.with_python
        self._cmake.definitions["BUILD_DOCUMENTATION"] = self.options.with_documentation
        self._cmake.configure()
        return self._cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy("LICENSE", dst="licenses")
        cmake = self._configure_cmake()
        cmake.install()
        tools.rmdir(os.path.join(self.package_folder, "lib", "cmake"))

    def package_info(self):
        self.cpp_info.names["cmake_find_package"] = "PROPOSAL"
        self.cpp_info.libs = ["PROPOSAL"]
        self.cpp_info.requires = [
            "cubicinterpolation::CubicInterpolation",
            "spdlog::spdlog",
            "nlohmann_json::nlohmann_json",
        ]

from conans import ConanFile, CMake


class PROPOSALConan(ConanFile):
    name = "PROPOSAL"
    version = "7.0.0"
    license = "GNU"
    author = "Maximilian Sackel mail@maxsac.de"
    url = "https://github.com/tudo-astroparticlephysics/PROPOSAL"
    description = "PROPOSAL (Propagator with optimal precision and optimized speed for all leptons) is presented as a public tool for propagating leptons and gamma rays through media. Up-to-date cross sections for ionization, bremsstrahlung, photonuclear interactions, electron pair production, Landau–Pomeranchuk–Migdal and Ter-Mikaelian effects, muon and tau decay, as well as Molière scattering are implemented for different parametrizations.  The full Paper can be found here."
    topics = ("lepton", "gamma", "propagator", "crosssection", "proposal")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }
    default_options = {"shared": False, "fPIC": True}
    generators = "cmake_find_package"
    exports_sources = "*"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def requirements(self):
        self.requires("CubicInterpolation/0.1")
        self.requires("spdlog/1.8.2")

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        self.copy("license*", dst="licenses", ignore_case=True, keep_path=False)
        self.copy("*.h", dst="include", src="src/PROPOSAL")
        self.copy("*.hpp", dst="include", src="src/PROPOSAL")
        self.copy("*.lib", dst="lib", keep_path=False)
        self.copy("*.dll", dst="bin", keep_path=False)
        self.copy("*.dylib*", dst="lib", keep_path=False)
        self.copy("*.so", dst="lib", keep_path=False)
        self.copy("*.a", dst="lib", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["PROPOSAL"]

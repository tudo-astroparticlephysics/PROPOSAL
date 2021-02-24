from conans import ConanFile, CMake, tools


class PROPOSALTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"

    def requirements(self):
        self.requires("gtest/1.10.0")

    def build(self):
        self.cmake = CMake(self)
        self.cmake.configure()
        self.cmake.build()

    def imports(self):
        self.copy("*.dll", dst="bin", src="bin")
        self.copy("*.dylib*", dst="bin", src="lib")
        self.copy("*.so*", dst="bin", src="lib")

    def test(self):
        self.run("ctest -j2 -E \(Brems\|Photo\|Epair\|Mupair\)")

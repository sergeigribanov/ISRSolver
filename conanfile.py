from conan import ConanFile
from conan.tools.files import copy
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps


class GeneratorRecipe(ConanFile):
    name = "isrsolver"
    version = "1.0"
    # Optional metadata
    license = "GPL-3.0"
    author = "Sergei Gribanov"
    url = "https://github.com/sergeigribanov/ISRSolver"
    description = "ISR tools"
    topics = ("ISR", "ISRSolver", "Kuraev", "Fading")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "*"

    def requirements(self):
        self.requires("nlohmann_json/3.11.2")
        self.requires("nlopt/2.7.1")
    
    def deploy(self):
        copy(self, "include/*", src=self.package_folder, dst=self.deploy_folder)
        copy(self, "lib/*", src=self.package_folder, dst=self.deploy_folder)
        copy(self, "bin/*", src=self.package_folder, dst=self.deploy_folder)
    
    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["BOOST_ROOT"] = "/sl/boost/Boost1.77.0_cc8_64/lib/cmake" # !!!
        tc.variables["GSL_INCLUDE_DIR"] = "/sl/gsl/gsl27_cc8_64/include" # !!!
        tc.variables["GSL_LIBRARY"] = "/sl/gsl/gsl27_cc8_64/lib/libgsl.so" # !!!
        tc.variables["GSL_CBLAS_LIBRARY"] =  "/sl/gsl/gsl27_cc8_64/lib/libgslcblas.so" # !!!
        tc.generate()
        
    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["isrsolver"]

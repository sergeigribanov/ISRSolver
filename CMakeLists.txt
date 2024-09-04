cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ISRSolver)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
find_package(Eigen3 REQUIRED NO_MODULE)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(NLOPT REQUIRED)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core MathCore Minuit2)
include(${ROOT_USE_FILE})
# find_package(Minuit2 REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package (Python3 COMPONENTS Development NumPy)
include_directories(${GSL_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${NLOPT_INCLUDE_DIRS})
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB PYSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/Py*.cpp)
list(FILTER SOURCES EXCLUDE REGEX "isrsolver.+.cpp")
list(FILTER SOURCES EXCLUDE REGEX "Py.+.cpp")
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
list(FILTER HEADERS EXCLUDE REGEX "Py.+.hpp")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/env.sh.in env.sh)

add_library(ISR SHARED ${SOURCES})
target_include_directories(ISR
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(ISR PUBLIC
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json)

add_executable(isrsolver-SLE ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLE.cpp)
target_link_libraries(isrsolver-SLE ISR)

add_executable(isrsolver-SLE2 ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLE2.cpp)
target_link_libraries(isrsolver-SLE2 ISR)

add_executable(isrsolver-calc-eff1d-2 ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-calc-eff1d-2.cpp)
target_link_libraries(isrsolver-calc-eff1d-2 ISR)

add_executable(isrsolver-Tikhonov ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov.cpp)
target_link_libraries(isrsolver-Tikhonov ISR)

add_executable(isrsolver-TSVD ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-TSVD.cpp)
target_link_libraries(isrsolver-TSVD ISR)

add_executable(isrsolver-calc-chi2
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-calc-chi2.cpp)
target_link_libraries(isrsolver-calc-chi2 ISR)

add_executable(isrsolver-SLE-chi2-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLE-chi2-test.cpp)
target_link_libraries(isrsolver-SLE-chi2-test ISR)

add_executable(isrsolver-iterative-chi2-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-iterative-chi2-test.cpp)
target_link_libraries(isrsolver-iterative-chi2-test ISR)

add_executable(isrsolver-Tikhonov-chi2-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-chi2-test.cpp)
target_link_libraries(isrsolver-Tikhonov-chi2-test ISR)

add_executable(isrsolver-SLE-ratio-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLE-ratio-test.cpp)
target_link_libraries(isrsolver-SLE-ratio-test ISR)

add_executable(isrsolver-Tikhonov-ratio-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-ratio-test.cpp)
target_link_libraries(isrsolver-Tikhonov-ratio-test ISR)

add_executable(isrsolver-condnum-test ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-condnum-test.cpp)
target_link_libraries(isrsolver-condnum-test ISR)

add_executable(isrsolver-draw-kernel ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-draw-kernel.cpp)
target_link_libraries(isrsolver-draw-kernel ISR)

add_executable(isrsolver-draw-kernel-integral ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-draw-kernel-integral.cpp)
target_link_libraries(isrsolver-draw-kernel-integral ISR)

add_executable(isrsolver-interp ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-interp.cpp)
target_link_libraries(isrsolver-interp ISR)

add_executable(isrsolver-LCurve-plot ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-LCurve-plot.cpp)
target_link_libraries(isrsolver-LCurve-plot ISR)

add_executable(isrsolver-Tikhonov-LCurve ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-LCurve.cpp)
target_link_libraries(isrsolver-Tikhonov-LCurve ISR)

add_executable(isrsolver-radcorr ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-radcorr.cpp)
target_link_libraries(isrsolver-radcorr ISR)

add_executable(isrsolver-radcorr-function ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-radcorr-function.cpp)
target_link_libraries(isrsolver-radcorr-function ISR)

add_executable(isrsolver-iterative ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-iterative.cpp)
target_link_libraries(isrsolver-iterative ISR)

add_executable(isrsolver-test-fit ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-test-fit.cpp)
target_link_libraries(isrsolver-test-fit ISR)

add_executable(isrsolver-iterative-use-vcs-fit-fcn ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-iterative-use-vcs-fit-fcn.cpp)
target_link_libraries(isrsolver-iterative-use-vcs-fit-fcn ISR)

add_executable(isrsolver-KuraevFadin-convolution
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-KuraevFadin-convolution.cpp)
target_link_libraries(isrsolver-KuraevFadin-convolution ISR)

add_library(PyISR SHARED ${PYSOURCES})
target_include_directories(PyISR
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>
  PRIVATE
  ${Python3_INCLUDE_DIRS}
  ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(PyISR ISR ${Python3_LIBRARIES})
set_target_properties(PyISR PROPERTIES PREFIX "")

set_target_properties(ISR
PROPERTIES
PUBLIC_HEADER "${HEADERS}")

install(TARGETS isrsolver-SLE DESTINATION bin)
install(TARGETS isrsolver-SLE2 DESTINATION bin)
install(TARGETS isrsolver-calc-eff1d-2 DESTINATION bin)
install(TARGETS isrsolver-Tikhonov DESTINATION bin)
install(TARGETS isrsolver-TSVD DESTINATION bin)
install(TARGETS isrsolver-LCurve-plot DESTINATION bin)
install(TARGETS isrsolver-Tikhonov-LCurve DESTINATION bin)
install(TARGETS isrsolver-calc-chi2 DESTINATION bin)
install(TARGETS isrsolver-SLE-chi2-test DESTINATION bin)
install(TARGETS isrsolver-iterative-chi2-test DESTINATION bin)
install(TARGETS isrsolver-Tikhonov-chi2-test DESTINATION bin)
install(TARGETS isrsolver-SLE-ratio-test DESTINATION bin)
install(TARGETS isrsolver-Tikhonov-ratio-test DESTINATION bin)
install(TARGETS isrsolver-condnum-test DESTINATION bin)
install(TARGETS isrsolver-iterative DESTINATION bin)
install(TARGETS isrsolver-iterative-use-vcs-fit-fcn DESTINATION bin)
install(TARGETS isrsolver-radcorr DESTINATION bin)
install(TARGETS isrsolver-radcorr-function DESTINATION bin)
install(TARGETS isrsolver-draw-kernel DESTINATION bin)
install(TARGETS isrsolver-draw-kernel-integral DESTINATION bin)
install(TARGETS isrsolver-KuraevFadin-convolution DESTINATION bin)
install(FILES ${PROJECT_BINARY_DIR}/env.sh DESTINATION bin)
install(TARGETS ISR
  EXPORT ISRSolverTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/ISRSolver)
install(TARGETS PyISR
  LIBRARY DESTINATION lib/python)

install(EXPORT ISRSolverTargets
  FILE ISRSolverTargets.cmake
  DESTINATION lib/cmake)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/ISRSolverConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ISRSolverConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/example"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ISRSolverConfigVersion.cmake"
  VERSION "${ISRSolver_VERSION_MAJOR}.${ISRSolver_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/ISRSolverConfig.cmake
  DESTINATION lib/cmake)
export(EXPORT ISRSolverTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/ISRSolverTargets.cmake")
export(PACKAGE ISRSolver)

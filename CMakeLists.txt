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
find_package(ROOT REQUIRED COMPONENTS RIO Net Core MathCore)
include(${ROOT_USE_FILE})
find_package(GSL REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(nlohmann_json REQUIRED)
include_directories(${GSL_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${NLOPT_INCLUDE_DIRS})
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLAE.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-TSVD.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-LCurve-plot.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-LCurve.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-calc-chi2.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLAE-chi2-test.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-chi2-test.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLAE-ratio-test.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-ratio-test.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-condnum-test.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-draw-kernel-integral.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-draw-kernel.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-interp.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-KuraevFadin-convolution.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-radcorr.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-iterative.cpp)
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/env.sh.in env.sh)

add_library(isr SHARED ${SOURCES})
target_include_directories(isr
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isr
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json)

add_executable(isrsolver-SLAE ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLAE.cpp)
target_include_directories(isrsolver-SLAE
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-SLAE
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-Tikhonov ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov.cpp)
target_include_directories(isrsolver-Tikhonov
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-Tikhonov
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-TSVD ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-TSVD.cpp)
target_include_directories(isrsolver-TSVD
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-TSVD
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-calc-chi2
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-calc-chi2.cpp)
target_include_directories(isrsolver-calc-chi2
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-calc-chi2
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-SLAE-chi2-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLAE-chi2-test.cpp)
target_include_directories(isrsolver-SLAE-chi2-test
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-SLAE-chi2-test
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-Tikhonov-chi2-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-chi2-test.cpp)
target_include_directories(isrsolver-Tikhonov-chi2-test
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-Tikhonov-chi2-test
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-SLAE-ratio-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-SLAE-ratio-test.cpp)
target_include_directories(isrsolver-SLAE-ratio-test
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-SLAE-ratio-test
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-Tikhonov-ratio-test
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-ratio-test.cpp)
target_include_directories(isrsolver-Tikhonov-ratio-test
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-Tikhonov-ratio-test
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-condnum-test ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-condnum-test.cpp)
target_include_directories(isrsolver-condnum-test
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-condnum-test
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-draw-kernel ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-draw-kernel.cpp)
target_include_directories(isrsolver-draw-kernel
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-draw-kernel
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(isrsolver-draw-kernel-integral ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-draw-kernel-integral.cpp)
target_include_directories(isrsolver-draw-kernel-integral
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-draw-kernel-integral
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-interp ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-interp.cpp)
target_include_directories(isrsolver-interp
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-interp
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-LCurve-plot ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-LCurve-plot.cpp)
target_include_directories(isrsolver-LCurve-plot
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-LCurve-plot
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-Tikhonov-LCurve ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-Tikhonov-LCurve.cpp)
target_include_directories(isrsolver-Tikhonov-LCurve
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-Tikhonov-LCurve
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-radcorr ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-radcorr.cpp)
target_include_directories(isrsolver-radcorr
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-radcorr
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(isrsolver-iterative ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-iterative.cpp)
target_include_directories(isrsolver-iterative
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-iterative
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)

add_executable(isrsolver-KuraevFadin-convolution
  ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-KuraevFadin-convolution.cpp)
target_include_directories(isrsolver-KuraevFadin-convolution
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-KuraevFadin-convolution
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
set_target_properties(isr
PROPERTIES
PUBLIC_HEADER "${HEADERS}")
install(TARGETS isrsolver-KuraevFadin-convolution DESTINATION bin)
install(TARGETS isrsolver-SLAE DESTINATION bin)
install(TARGETS isrsolver-TSVD DESTINATION bin)
install(TARGETS isrsolver-Tikhonov DESTINATION bin)
install(TARGETS isrsolver-condnum-test DESTINATION bin)
install(TARGETS isrsolver-iterative DESTINATION bin)
install(TARGETS isrsolver-Tikhonov-LCurve DESTINATION bin)
install(TARGETS isrsolver-radcorr DESTINATION bin)
install(TARGETS isrsolver-LCurve-plot DESTINATION bin)
install(TARGETS isrsolver-draw-kernel DESTINATION bin)
install(TARGETS isrsolver-draw-kernel-integral DESTINATION bin)
install(FILES ${PROJECT_BINARY_DIR}/env.sh DESTINATION bin)
install(TARGETS isr
  EXPORT ISRSolverTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/ISRSolver)
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

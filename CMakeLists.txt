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
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-bias-test.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-chi2-test.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/test-draw-kernel-integral.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/test-draw-kernel.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/test-interp.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/test-chi2.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/convolution.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/lcurve.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-lcurve.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/reg-perturb-errors.cpp)
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
add_executable(isrsolver ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_include_directories(isrsolver
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(isrsolver-chi2-test ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-chi2-test.cpp)
target_include_directories(isrsolver-chi2-test
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-chi2-test
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(isrsolver-bias-test ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-bias-test.cpp)
target_include_directories(isrsolver-bias-test
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-bias-test
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(test-draw-kernel ${CMAKE_CURRENT_SOURCE_DIR}/src/test-draw-kernel.cpp)
target_include_directories(test-draw-kernel
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(test-draw-kernel
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(test-draw-kernel-integral ${CMAKE_CURRENT_SOURCE_DIR}/src/test-draw-kernel-integral.cpp)
target_include_directories(test-draw-kernel-integral
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(test-draw-kernel-integral
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(test-interp ${CMAKE_CURRENT_SOURCE_DIR}/src/test-interp.cpp)
target_include_directories(test-interp
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(test-interp
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(test-chi2 ${CMAKE_CURRENT_SOURCE_DIR}/src/test-chi2.cpp)
target_include_directories(test-chi2
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(test-chi2
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)



add_executable(lcurve ${CMAKE_CURRENT_SOURCE_DIR}/src/lcurve.cpp)
target_include_directories(lcurve
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(lcurve
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(isrsolver-lcurve ${CMAKE_CURRENT_SOURCE_DIR}/src/isrsolver-lcurve.cpp)
target_include_directories(isrsolver-lcurve
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrsolver-lcurve
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(reg-perturb-errors ${CMAKE_CURRENT_SOURCE_DIR}/src/reg-perturb-errors.cpp)
target_include_directories(reg-perturb-errors
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(reg-perturb-errors
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${ROOT_LIBRARIES}
  nlohmann_json::nlohmann_json
  isr)
add_executable(isrconv ${CMAKE_CURRENT_SOURCE_DIR}/src/convolution.cpp)
target_include_directories(isrconv
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/ISRSolver>)
target_link_libraries(isrconv
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
install(TARGETS isrsolver DESTINATION bin)
install(TARGETS lcurve DESTINATION bin)
install(TARGETS isrconv DESTINATION bin)
install(TARGETS isrsolver-lcurve DESTINATION bin)
install(TARGETS reg-perturb-errors DESTINATION bin)
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

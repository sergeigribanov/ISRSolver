cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ISRSolver)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core)
include(${ROOT_USE_FILE})
find_package(GSL REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(isrsolver ${SOURCES})
target_include_directories(isrsolver
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  ${GSL_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS})
target_link_libraries(isrsolver
  PRIVATE
  ${ROOT_LIBRARIES}
  GSL::gsl GSL::gslcblas
  ${Boost_LIBRARIES})
install(TARGETS isrsolver
  RUNTIME DESTINATION bin)

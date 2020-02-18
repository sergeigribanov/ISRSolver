cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ISRSolver)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
find_package(Eigen3 REQUIRED NO_MODULE)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core)
include(${ROOT_USE_FILE})
find_package(GSL REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/convolution.cpp)
file(GLOB CONV_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM CONV_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
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
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${ROOT_LIBRARIES})
add_executable(isrconv ${CONV_SOURCES})
target_include_directories(isrconv
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  ${GSL_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS})
target_link_libraries(isrconv
  PRIVATE
  GSL::gsl GSL::gslcblas
  Eigen3::Eigen
  ${Boost_LIBRARIES}
  ${ROOT_LIBRARIES})
install(TARGETS isrsolver DESTINATION bin)
install(TARGETS isrconv DESTINATION bin)

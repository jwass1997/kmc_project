cmake_minimum_required(VERSION 3.10)
project(kmc_project VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_ROOT ${CMAKE_SOURCE_DIR}) # .../kmc_project
set(EXTERNAL_ROOT ${CMAKE_SOURCE_DIR}/../local) #parent dir: ../local

set(BOOST_INCLUDEDIR $ENV{BOOST_INC_DIR})
set(BOOST_LIBRARYDIR $ENV{BOOST_LIB_DIR})

find_package(Boost REQUIRED COMPONENTS program_options)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find Boost. Did you set BOOST_LIB_DIR and BOOST_LIB_DIR?")
endif()

find_package(ZLIB REQUIRED)
if(NOT ZLIB_FOUND)
  message(FATAL_ERROR "Could not find zlib on your system")
endif()

find_package(OpenMP REQUIRED)
if(NOT OpenMP_CXX_FOUND)
  message(FATAL_ERROR "Could not find an OpenMP-capable compiler")
endif()

find_library(CNPY_LIB NAMES cnpy PATHS ${EXTERNAL_ROOT}/lib NO_DEFAULT_PATH)
if(NOT CNPY_LIB)
  message(FATAL_ERROR
    "Could not find cnpy library in ${EXTERNAL_ROOT}/lib")
endif()

find_library(MFEM_LIB NAMES mfem PATHS ${EXTERNAL_ROOT}/lib NO_DEFAULT_PATH)
if(NOT MFEM_LIB)
  message(FATAL_ERROR
    "Could not find mfem library in ${EXTERNAL_ROOT}/lib")
endif()

set(SOURCES
    src/main.cpp
    src/State.cpp 
    src/KMCSimulator.cpp 
    src/Configuration.cpp
    src/utils.cpp 
    src/FEMmethods.cpp 
    src/Random.cpp 
)

add_executable(kmc_project ${SOURCES})

target_include_directories(kmc_project PRIVATE
    ${PROJECT_ROOT}/include
    ${EXTERNAL_ROOT}/include
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS} 
)

target_compile_options(kmc_project PRIVATE 
    -Ofast
    ${OpenMP_CXX_FLAGS}
)

target_link_libraries(kmc_project PRIVATE
    ${CNPY_LIB}
    ${MFEM_LIB}
    ${Boost_LIBRARIES}
    ${ZLIB_LIBRARIES}
    OpenMP::OpenMP_CXX
)
# basic setup for cmake
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD_REQUIRED 17)
set(CMAKE_CXX_STANDARD 17)

set(DEFAULT_BUILD "RELEASE")
set(VALID_BUILDS DEVEL RELEASE DEBUG)
set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD}" CACHE STRING  "Choose the type of build, options are: ${VALID_BUILDS}")
string(TOUPPER "${CMAKE_BUILD_TYPE}" CHECK_BUILD)
LIST(FIND VALID_BUILDS "${CHECK_BUILD}" CURRENT_BUILD)

if(CURRENT_BUILD LESS 0)
  message(SEND_ERROR "Please specify a valid CMAKE_BUILD_TYPE type! Valid build types are: " "${VALID_BUILDS}")
endif()

set(CMAKE_CUDA_ARCHITECTURES "native")

#
# PROJECT is Axiton
#
project("AXITON")

#######################################################################
# Axions OPTIONS
# likely to be changed by users
#######################################################################

set(CMAKE_CUDA "/usr/local/cuda")
#set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
set(CMAKE_CUDA_COMPILER "$ENV{CUDA_HOME}/bin/nvcc")
enable_language(CUDA)
# do all the build definitions

include_directories(SYSTEM ${HDF5_INCLUDE_DIR})

find_package(OpenMP REQUIRED)
find_package(HDF5   REQUIRED)

#find_library(FFTW_LIBRARY      NAMES fftw3 fftw           HINTS "$ENV{FFTW_HOME}")
#find_library(FFTW_OMP_LIBRARY  NAMES fftw3_omp fftw_omp   PATHS $ENV{FFTW_HOME}/lib)

set(FFTW_HOME $ENV{FFTW_HOME})
find_library(FFTW_LIBRARY     NAMES fftw3 fftw         HINTS ${FFTW_HOME} PATH_SUFFIXES lib lib64)

# Detect compiler for selecting the correct OpenMP-enabled FFTW library
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(_fftw_omp_name fftw3_omp_GNU)
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(_fftw_omp_name fftw3_omp_Intel)
elseif(CMAKE_C_COMPILER_ID MATCHES "NVIDIA")
    set(_fftw_omp_name fftw3_omp_nvidia)
else()
    set(_fftw_omp_name fftw3_omp) # Generic case
endif()

find_library(FFTW_OMP_LIBRARY NAMES ${_fftw_omp_name} HINTS ${FFTW_HOME} PATH_SUFFIXES lib lib64)

if(NOT FFTW_OMP_LIBRARY)
message(FATAL_ERROR "Could not find OpenMP FFTW library for compiler '${_fftw_omp_name}' in ${FFTW_HOME}")
endif()


set(GPU_ARCH sm_70 CACHE STRING "set the GPU architecture (sm_50, sm_52, sm_60, sm_61, sm_62, sm_70)")
set(GPU_COMP compute_70 CACHE STRING "set the GPU compute capability (compute_50, compute_52, compute_60, compute_61, compute_62, compute_70)")
set(DEFARCH "700")
#find_package(CUDA REQUIRED)
#set(CUDA_HOST_COMPILER gcc)
set(CUDA_PROPAGATE_HOST_FLAGS ON)
set(CUDA_NVCC_FLAGS "-arch=${GPU_COMP} -code=${GPU_ARCH} -ftz=true -prec-div=false -prec-sqrt=false -std=c++17 --expt-relaxed-constexpr")

set(CUDA_NVCC_FLAGS_DEVEL "-O3 -lineinfo" CACHE STRING
    "Flags used by the CUDA compiler during regular development builds." FORCE )
set(CUDA_NVCC_FLAGS_RELEASE "-O3 -w" CACHE STRING
  "Flags used by the C++ compiler during release builds." FORCE )
set(CUDA_NVCC_FLAGS_HOSTDEBUG "-g -G" CACHE STRING
    "Flags used by the C++ compiler during debug builds." FORCE )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -mcmodel=medium -std=c++17 -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -mcmodel=medium -std=c++17 -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wa,-q -lc++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q -lc++")
endif()

set(CMAKE_CXX_FLAGS_DEVEL "-O3 -Wall" CACHE STRING
    "Flags used by the C++ compiler during regular development builds." FORCE )
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -w" CACHE STRING
    "Flags used by the C++ compiler during release builds." FORCE )
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -fno-inline" CACHE STRING
    "Flags used by the C++ compiler during debug builds." FORCE )

include_directories("$ENV{HDF5_HOME}/include")
include_directories("$ENV{FFTW_HOME}/include")
include_directories(include)
include_directories(SYSTEM ${CUDA_INCLUDE_DIRS})

### add tests and library
add_subdirectory(lib)
add_subdirectory(src)

cmake_minimum_required(VERSION 3.8)
project(nebula CXX CUDA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# Find HDF5 library
find_package(HDF5 1.8.12 REQUIRED CXX HL)
include_directories(${HDF5_INCLUDE_DIRS})
list(APPEND LIBRARIES ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})

# Find MPI library if compiling against parallel HDF5
if(HDF5_IS_PARALLEL)
	find_package(MPI REQUIRED)
	include_directories(${MPI_CXX_INCLUDE_PATH})
	list(APPEND LIBRARIES ${MPI_CXX_LIBRARIES})
endif()

# Find threading library
find_package(Threads REQUIRED)


string(APPEND CMAKE_CXX_FLAGS " -O3")
string(APPEND CMAKE_CUDA_FLAGS " --expt-relaxed-constexpr")
string(APPEND CMAKE_CUDA_FLAGS " -O3 -prec-sqrt=false")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_30,code=compute_30")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_37,code=compute_37")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_60,code=compute_60")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_61,code=compute_61")
if(CMAKE_CUDA_COMPILER_VERSION LESS 9.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_20,code=compute_20")
endif(CMAKE_CUDA_COMPILER_VERSION LESS 9.0)

include_directories(3rdparty)
add_subdirectory(source)

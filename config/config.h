#ifndef __CONFIG_H_
#define __CONFIG_H_

// Settings
#define USE_DOUBLE 0

// Include relevant headers
#include <cstdint>
#include <cstddef>
#include <stdexcept>

#include "device.h"
#include "data_types.h"

#include "../common/scalar_math.h"
#include "../common/vector_math.h"

#if CUDA_AVAILABLE
	#include "../common/cuda/cuda_new.h"
	#include "../common/cuda/cuda_mem_scope.h"
#endif // CUDA_AVAILABLE

// TODO: Put these somewhere nicer.
// Define range and resolution for cross section tables.
GLOBAL_CONSTANT real K_min = 1;    // eV
GLOBAL_CONSTANT real K_max = 50e3; // eV
GLOBAL_CONSTANT int K_cnt = 1024;
GLOBAL_CONSTANT int P_cnt = 1024;

// A constant of mathematics
GLOBAL_CONSTANT real pi = 3.1415926535897932_r;

#endif // __CONFIG_H_

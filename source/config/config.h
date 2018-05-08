#ifndef __CONFIG_H_
#define __CONFIG_H_

/**
 * \file config/config.h
 */

/**
 * \brief Whether or not to use double precision in simulations.
 *
 * Set to 0 for single-precision and to 1 for double-precision calculations.
 */
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
constexpr real K_min = 1;    ///< Min energy in simulation, in eV
constexpr real K_max = 50e3; ///< Max energy in simulation, in eV
constexpr int K_cnt = 1024;  ///< Number of data points along energy axis
constexpr int P_cnt = 1024;  ///< Number of data points along probability axis

/// A frequently-used constant of mathematics
constexpr real pi = 3.1415926535897932_r;

#endif // __CONFIG_H_

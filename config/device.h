#ifndef __DEVICE_H_
#define __DEVICE_H_

/**
 * \file config/device.h
 * \brief Global definitions for GPU code
 */

/**
 * \def CUDA_AVAILABLE
 * \brief Set to 1 if CUDA is available, and to zero otherwise.
 *
 * This is not a setting! We just check if the CUDA compiler is being used.
 */
#ifdef __CUDACC__
	#define CUDA_AVAILABLE 1
#else
	#define CUDA_AVAILABLE 0
#endif

/**
 * \def CUDA_COMPILING
 * \brief Detect if we are compiling for CUDA right now.
 *
 * This is not a setting! It's just a check if CUDA is the target right now.
 */
#ifdef __CUDA_ARCH__
	#define CUDA_COMPILING 1
#else
	#define CUDA_COMPILING 0
#endif

/**
 * \def CPU
 * \brief Designates that a given function must be run on a CPU, and never on a GPU.
 */

/**
 * \def GPU
 * \brief Designates that a given function must be run on a GPU, and never on a CPU.
 */

/**
 * \def PHYSICS
 * \brief Designates that a given function may either be run on a CPU or a GPU.
 * This is common for the physics code, hence the name.
 */
#if CUDA_AVAILABLE
	// Load CUDA headers
	#include <cuda_runtime.h>
	#include <device_launch_parameters.h>
	#include <device_functions.h>

	#define CPU __host__                // CPU-only code
	#define GPU __device__              // GPU-only code
	#define PHYSICS __host__ __device__ // Physics are GPU and CPU code
#else // CUDA_AVAILABLE
	#define CPU
	#define GPU
	#define PHYSICS
#endif

#endif // __DEVICE_H_

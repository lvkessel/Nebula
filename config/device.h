#ifndef __DEVICE_H_
#define __DEVICE_H_

/*
 * Global definitions for GPU code
 */

// Detect availability of CUDA functions
#ifdef __CUDACC__
	#define CUDA_AVAILABLE 1
#else
	#define CUDA_AVAILABLE 0
#endif

// Detect if we are compiling for CUDA right now
#ifdef __CUDA_ARCH__
	#define CUDA_COMPILING 1
#else
	#define CUDA_COMPILING 0
#endif

#if CUDA_AVAILABLE
	// Load CUDA headers
	#include <cuda_runtime.h>
	#include <device_launch_parameters.h>
	#include <device_functions.h>

	#define HOST __host__               // CPU-only code
	#define PHYSICS __host__ __device__ // Physics are GPU and CPU code
#else // CUDA_AVAILABLE
	#define HOST
	#define PHYSICS
#endif

#if CUDA_COMPILING
	#define GLOBAL_CONSTANT const __device__
#else
	#define GLOBAL_CONSTANT const
#endif

#endif // __DEVICE_H_
#ifndef __DATA_TYPES_H_
#define __DATA_TYPES_H_

/*
 * Define data types to be used in device code
 */

// Index type
using std::size_t;

// scalar and vector types
#if CUDA_AVAILABLE
	#if USE_DOUBLE
		using real = double;
		using vec3 = double3;
		using vec4 = double4;
	#else // USE_DOUBLE
		using real = float;
		using vec3 = float3;
		using vec4 = float4;
	#endif // USE_DOUBLE
#else // CUDA_AVAILABLE
	#if USE_DOUBLE
		using real = double;
		struct vec3 { double x, y, z; };
		struct vec4 { double x, y, z, w; };
	#else // USE_DOUBLE
		using real = float;
		struct vec3 { float x, y, z; };
		struct vec4 { float x, y, z, w; };
	#endif // USE_DOUBLE
#endif // CUDA_AVAILABLE

// Small value limiting accuracy of numerical computations
#include <limits>
GLOBAL_CONSTANT real EPSILON = 10 * std::numeric_limits<real>::epsilon();

// User-defined literal for using reals: e.g. 1.32_r == (real)1.32.
constexpr PHYSICS real operator"" _r(long double r)
{
	return static_cast<real>(r);
}

#endif // __DATA_TYPES_H_
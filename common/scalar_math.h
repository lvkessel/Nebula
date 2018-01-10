#ifndef __SCALAR_MATH_H_
#define __SCALAR_MATH_H_

/*
 * Scalar math to be used in device code
 */

inline PHYSICS int clampi(int i, int low, int high) { return i < low ? low : (i > high ? high : i); }
inline PHYSICS real clampr(real i, real low, real high) { return i < low ? low : (i > high ? high : i); }

#if CUDA_COMPILING
	#include <math_functions.h>
	#include <math_constants.h>

	#if USE_DOUBLE
		inline __device__ int rintr(real i) { return __double2int_rn(i); }
		inline __device__ int floorr(real i) { return __double2int_rd(i); }
		inline __device__ real minr(real a, real b) { return fmin(a, b); }
		inline __device__ real maxr(real a, real b) { return fmax(a, b); }
		inline __device__ real absr(real x) { return fabs(x); }
		inline __device__ real sqrtr(real x) { return sqrt(x); }
		inline __device__ real cbrtr(real x) { return cbrt(x); }
		inline __device__ real logr(real x) { return log(x); }
		inline __device__ real expr(real x) { return exp(x); }
		inline __device__ real expm1r(real x) { return expm1(x); }
		inline __device__ real sinr(real x) { return sin(x); }
		inline __device__ real cosr(real x) { return cos(x); }
		inline __device__ real atan2r(real y, real x) { return atan2(y, x); }
		inline __device__ void sincosr(real x, real* sin_ptr, real* cos_ptr) { sincos(x, sin_ptr, cos_ptr); }
		inline __device__ real copysignr(real x, real y) { return copysign(x, y); }
		inline __device__ real saturater(real x) { return clampr(x, 0, 1); }
	#else // USE_DOUBLE
		inline __device__ int rintr(real i) { return __float2int_rn(i); }
		inline __device__ int floorr(real i) { return __float2int_rd(i); }
		inline __device__ real minr(real a, real b) { return fminf(a, b); }
		inline __device__ real maxr(real a, real b) { return fmaxf(a, b); }
		inline __device__ real absr(real x) { return fabsf(x); }
		inline __device__ real sqrtr(real x) { return sqrtf(x); }
		inline __device__ real cbrtr(real x) { return cbrtf(x); }
		inline __device__ real logr(real x) { return logf(x); }
		inline __device__ real expr(real x) { return expf(x); }
		inline __device__ real expm1r(real x) { return expm1f(x); }
		inline __device__ real sinr(real x) { return sinf(x); }
		inline __device__ real cosr(real x) { return cosf(x); }
		inline __device__ real atan2r(real y, real x) { return atan2f(y, x); }
		inline __device__ void sincosr(real x, real* sin_ptr, real* cos_ptr) { sincosf(x, sin_ptr, cos_ptr); }
		inline __device__ real copysignr(real x, real y) { return copysignf(x, y); }
		inline __device__ real saturater(real x) { return __saturatef(x); }
	#endif // USE_DOUBLE

#else // CUDA_COMPILING
	#include <cmath>
	#include <algorithm>

	inline CPU int rintr(real i) { return (int)std::rint(i); }
	inline CPU int floorr(real i) { return (int)std::floor(i); }
	inline CPU real minr(real a, real b) { return std::min(a, b); }
	inline CPU real maxr(real a, real b) { return std::max(a, b); }
	inline CPU real absr(real x) { return std::abs(x); }
	inline CPU real sqrtr(real x) { return std::sqrt(x); }
	inline CPU real cbrtr(real x) { return std::cbrt(x); }
	inline CPU real logr(real x) { return std::log(x); }
	inline CPU real expr(real x) { return std::exp(x); }
	inline CPU real expm1r(real x) { return std::expm1(x); }
	inline CPU real sinr(real x) { return std::sin(x); }
	inline CPU real cosr(real x) { return std::cos(x); }
	inline CPU real atan2r(real y, real x) { return std::atan2(y, x); }

	inline CPU void sincosr(real x, real* sin_ptr, real* cos_ptr)
	{
		// TODO: try using sin^2(x) + cos^2(x) = 1,
		// but remember to get the sign right.
		*sin_ptr = sinr(x);
		*cos_ptr = cosr(x);
	}
	inline CPU real copysignr(real x, real y) { return std::copysign(x, y); }
	inline CPU real saturater(real x) { return clampr(x, 0, 1); }
#endif // CUDA_COMPILING

#endif // __SCALAR_MATH_H_

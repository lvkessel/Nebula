#ifndef __SCALAR_MATH_H_
#define __SCALAR_MATH_H_

#include <cmath>
#include <algorithm>
#if CUDA_AVAILABLE
	#include <math_functions.h>
	#include <math_constants.h>
#endif // CUDA_AVAILABLE

/*
 * Scalar math to be used in device code
 */

inline PHYSICS int clampi(int i, int low, int high)
{
	return i < low ? low : (i > high ? high : i);
}
inline PHYSICS real clampr(real i, real low, real high)
{
	return i < low ? low : (i > high ? high : i);
}

inline PHYSICS int rintr(real i)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return __double2int_rn(i);
	#else // USE_DOUBLE
		return __float2int_rn(i);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return (int)std::rint(i);
#endif // CUDA_COMPILING
}

inline PHYSICS int floorr(real i)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return __double2int_rd(i);
	#else // USE_DOUBLE
		return __float2int_rd(i);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return (int)std::floor(i);
#endif // CUDA_COMPILING
}

inline PHYSICS real minr(real a, real b)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return fmin(a, b);
	#else // USE_DOUBLE
		return fminf(a, b);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::min(a, b);
#endif // CUDA_COMPILING
}

inline PHYSICS real maxr(real a, real b)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return fmax(a, b);
	#else // USE_DOUBLE
		return fmaxf(a, b);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::max(a, b);
#endif // CUDA_COMPILING
}

inline PHYSICS real absr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return fabs(x);
	#else // USE_DOUBLE
		return fabsf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::abs(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real sqrtr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return sqrt(x);
	#else // USE_DOUBLE
		return sqrtf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::sqrt(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real cbrtr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return cbrt(x);
	#else // USE_DOUBLE
		return cbrtf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::cbrt(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real logr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return log(x);
	#else // USE_DOUBLE
		return logf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::log(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real expr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return exp(x);
	#else // USE_DOUBLE
		return expf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::exp(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real expm1r(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return expm1(x);
	#else // USE_DOUBLE
		return expm1f(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::expm1(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real sinr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return sin(x);
	#else // USE_DOUBLE
		return sinf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::sin(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real cosr(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return cos(x);
	#else // USE_DOUBLE
		return cosf(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::cos(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real atan2r(real y, real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return atan2(y, x);
	#else // USE_DOUBLE
		return atan2f(y, x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	 return std::atan2(y, x);
#endif // CUDA_COMPILING
}

inline PHYSICS void sincosr(real x, real* sin_ptr, real* cos_ptr)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		sincos(x, sin_ptr, cos_ptr);
	#else // USE_DOUBLE
		sincosf(x, sin_ptr, cos_ptr);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	// TODO: try using sin^2(x) + cos^2(x) = 1,
	// but remember to get the sign right.
	*sin_ptr = std::sin(x);
	*cos_ptr = std::cos(x);
#endif // CUDA_COMPILING
}

inline PHYSICS real copysignr(real x, real y)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return copysign(x, y);
	#else // USE_DOUBLE
		return copysignf(x, y);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	return std::copysign(x, y);
#endif // CUDA_COMPILING
}

inline PHYSICS real saturater(real x)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return clampr(x, 0, 1);
	#else // USE_DOUBLE
		return __saturatef(x);
	#endif // USE_DOUBLE
#else // CUDA_COMPILING
	return clampr(x, 0, 1);
#endif // CUDA_COMPILING
}

#endif // __SCALAR_MATH_H_

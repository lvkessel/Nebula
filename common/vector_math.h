#ifndef __VECTOR_MATH_H_
#define __VECTOR_MATH_H_

/*
 * Vector math to be used in host or device code
 */

inline PHYSICS vec3 operator+(vec3 a, vec3 b)
{
	return{ a.x + b.x, a.y + b.y, a.z + b.z };
}

inline PHYSICS vec3 operator-(vec3 a, vec3 b)
{
	return{ a.x - b.x, a.y - b.y, a.z - b.z };
}

inline PHYSICS vec3 operator*(vec3 a, real b)
{
	return{ a.x*b, a.y*b, a.z*b };
}

inline PHYSICS vec3 operator*(real a, vec3 b)
{
	return{ a*b.x, a*b.y, a*b.z };
}

inline PHYSICS vec3 operator/(vec3 a, real b)
{
	return a * (1 / b);
}

inline PHYSICS vec3 operator-(vec3 a)
{
	return{ -a.x, -a.y, -a.z };
}

inline PHYSICS void operator+=(vec3& a, vec3 b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

inline PHYSICS void operator-=(vec3& a, vec3 b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

inline PHYSICS void operator*=(vec3& a, real b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
}

inline PHYSICS void operator/=(vec3& a, real b)
{
	a *= 1 / b;
}

inline PHYSICS vec3 cross_product(vec3 a, vec3 b)
{
	return
	{
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	};
}

inline PHYSICS real dot_product(vec3 a, vec3 b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline PHYSICS real magnitude(vec3 a)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		return norm3d(a.x, a.y, a.z);
	#else // USE_DOUBLE
		return norm3df(a.x, a.y, a.z);
	#endif //USE_DOUBLE
#else // CUDA_COMPILING
	return sqrtr(a.x*a.x + a.y*a.y + a.z*a.z);
#endif // CUDA_COMPILING
}

inline PHYSICS real magnitude_squared(vec3 a)
{
	return a.x*a.x + a.y*a.y + a.z*a.z;
}

inline PHYSICS vec3 normalised(vec3 a)
{
#if CUDA_COMPILING
	#if USE_DOUBLE
		const real rnorm = rnorm3d(a.x, a.y, a.z);
	#else // USE_DOUBLE
		const real rnorm = rnorm3df(a.x, a.y, a.z);
	#endif //USE_DOUBLE
	return a * rnorm;
#else // CUDA_COMPILING
	return a / sqrtr(a.x*a.x + a.y*a.y + a.z*a.z);
#endif // CUDA_COMPILING
}

/*
 * Find a vector normal to dir.
 * This function finds two axes normal to dir and returns a vector with angle phi in that system.
 * If the input direction is is a unit vector, the returned vector is also a unit vector.
 * Otherwise, the returned vector is not a unit either but it is normal to the input direction.
 */
inline PHYSICS vec3 make_normal_vec(vec3 dir, real phi)
{
	real sin_azimuth, cos_azimuth;
	sincosr(atan2r(dir.y, dir.x), &sin_azimuth, &cos_azimuth);

	const vec3 unit_v {
		dir.z*cos_azimuth,
		dir.z*sin_azimuth,
		-sqrtr(dir.x*dir.x + dir.y*dir.y)
	};
	const vec3 unit_u = cross_product(unit_v, dir);

	real sin_phi, cos_phi;
	sincosr(phi, &sin_phi, &cos_phi);
	return unit_u*cos_phi + unit_v*sin_phi;
}

/*
 * Create a unit vector, given Euler angles
 */
inline PHYSICS vec3 make_unit_vec(real cos_theta, real phi)
{
	real sin_phi, cos_phi;
	sincosr(phi, &sin_phi, &cos_phi);

	cos_theta = clampr(cos_theta, -1, 1);
	const real sin_theta = sqrtf(1 - cos_theta*cos_theta);
	return vec3 {
		sin_theta * cos_phi,
		sin_theta * sin_phi,
		cos_theta
	};
}

#endif // VECTOR_MATH_H_

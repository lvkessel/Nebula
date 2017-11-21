HOST triangle::triangle(vec3 r0, vec3 r1, vec3 r2, int material_in, int material_out)
	: material_in(material_in), material_out(material_out),
	_r0(r0), _e1(r1-r0), _e2(r2-r0)
{}

// TODO: trayce implementation
PHYSICS real triangle::intersect_ray(vec3 ray_start, vec3 ray_direction) const
{
	const vec3 pvec = cross_product(ray_direction, _e2);
	const real det = dot_product(_e1, pvec);
	if(absr(det) < EPSILON)
		return -1;

	const real u = dot_product(ray_start - _r0, pvec) / det;
	if((u < -EPSILON) || (u > 1 + EPSILON))
		return -1;

	const vec3 qvec = cross_product(ray_start - _r0, _e1);
	const real v = dot_product(ray_direction, qvec) / det;
	if((v < -EPSILON) || (u+v > 1+EPSILON))
		return -1;

	const real t = dot_product(_e2, qvec) / det;
	if(t < 0)
		return -1;

	return t;
}

PHYSICS vec3 triangle::r0() const
{
	return _r0;
}
PHYSICS vec3 triangle::r1() const
{
	return _r0 + _e1;
}
PHYSICS vec3 triangle::r2() const
{
	return _r0 + _e2;
}

PHYSICS vec3 triangle::AABB_min() const
{
	const vec3 A = r0();
	const vec3 B = r1();
	const vec3 C = r2();

	return vec3
	{
		minr(minr(A.x, B.x), C.x),
		minr(minr(A.y, B.y), C.y),
		minr(minr(A.z, B.z), C.z),
	};
}

PHYSICS vec3 triangle::AABB_max() const
{
	const vec3 A = r0();
	const vec3 B = r1();
	const vec3 C = r2();

	return vec3
	{
		maxr(maxr(A.x, B.x), C.x),
		maxr(maxr(A.y, B.y), C.y),
		maxr(maxr(A.z, B.z), C.z),
	};
}

PHYSICS vec3 triangle::get_normal() const
{
	return cross_product(_e1, _e2);
}
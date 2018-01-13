#ifndef __TRIANGLE_H_
#define __TRIANGLE_H_

class triangle
{
public:
	// Constructor, takes vertices for the edges.
	// TODO: document definition for "in" and "out" materials
	CPU triangle(vec3 r0, vec3 r1, vec3 r2, int material_in, int material_out);

	// Next ray-triangle intersection is at ray_start + t*ray_direction,
	// where t is the return value of this function.
	// Returns -1 if there is no intersection.
	PHYSICS real intersect_ray(vec3 ray_start, vec3 ray_direction) const;

	// Get vertices (r1 and r2 require computation)
	PHYSICS vec3 r0() const;
	PHYSICS vec3 r1() const;
	PHYSICS vec3 r2() const;

	// Get bounding box (requires computation)
	PHYSICS vec3 AABB_min() const;
	PHYSICS vec3 AABB_max() const;

	// Get unnormalised normal (requires computation)
	PHYSICS vec3 get_normal() const;

	// TODO: typename material_manager::material_index_t
	int material_in;
	int material_out;

private:
	vec3 _r0;
	vec3 _e1;
	vec3 _e2;
};

#include "triangle.inl"

#endif // __TRIANGLE_H_

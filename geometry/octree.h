#ifndef __GEOMETRY_OCTREE_H_
#define __GEOMETRY_OCTREE_H_

/*
 * Geometry class.
 * This class is responsible for the collision detection system.
 * It holds the simulation domain and the triangles.
 * 
 * This is an octree implementation.
 */

#include "../core/triangle.h"
#include "../core/events.h"

namespace nbl { namespace geometry {

template<bool gpu_flag>
struct _octree_factory;

template<bool gpu_flag>
class octree
{
public:
	using triangle_index_t = uint32_t;

	static HOST octree<gpu_flag> create(std::vector<triangle> const & triangles);
	static HOST void destroy(octree & geometry);

	// Returns whether a certain position is part of the simulation domain
	inline PHYSICS bool in_domain(vec3 position) const;

	// Try to propagate from start to start + distance*direction.
	// The triangle_event contains a nullptr triangle pointer if nothing is in the way.
	// An optional ignore_triangle may be passed to indicate a single triangle that must be ignored.
	// TODO ignore_material's datatype
	inline PHYSICS intersect_event propagate(vec3 start, vec3 direction, real distance,
		triangle const * ignore_triangle, int ignore_material);

	// Get the maximum distance that can be travelled inside the simulation domain.
	inline PHYSICS real get_max_extent() const;

	// Get the (axis-aligned) simulation domain
	inline PHYSICS vec3 AABB_min() const;
	inline PHYSICS vec3 AABB_max() const;

private:
	inline HOST void set_AABB(vec3 min, vec3 max);
	inline static PHYSICS vec3 AABB_intersect(vec3 pos, vec3 dir, vec3 center, vec3 halfsize);

	inline static PHYSICS int clz(uint64_t x);

	int* _octree_data    = nullptr;
	triangle* _triangles = nullptr;
	triangle_index_t _N  = 0;
	vec3 _AABB_center    = { 0, 0, 0 };
	vec3 _AABB_halfsize  = { 0, 0, 0 };
	real _max_extent     = 0;

	friend struct _octree_factory<gpu_flag>;
};

}} // namespace nbl::geometry

#include "octree.inl"

#endif // __GEOMETRY_OCTREE_H_
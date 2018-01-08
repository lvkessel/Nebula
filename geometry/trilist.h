#ifndef __GEOMETRY_TRILIST_H_
#define __GEOMETRY_TRILIST_H_

/*
 * Geometry class.
 * This class is responsible for the collision detection system.
 * It holds the simulation domain and the triangles.
 * 
 * This is the simplest implementation, which simply holds a list of
 * all triangles and tries all of them in turn when checking for collisions.
 */

#include "../core/triangle.h"
#include "../core/events.h"

namespace nbl { namespace geometry {

template<bool gpu_flag>
struct _trilist_factory;

template<bool gpu_flag>
class trilist
{
public:
	using triangle_index_t = uint32_t;

	static CPU trilist create(std::vector<triangle> const & triangles);
	static CPU void destroy(trilist & geometry);

	// Returns whether a certain position is part of the simulation domain
	inline PHYSICS bool in_domain(vec3 position);

	// Try to propagate from start to start+direction.
	// The triangle_event contains a nullptr triangle pointer if nothing is in the way.
	// An optional ignore_triangle may be passed to indicate a single triangle that must be ignored.
	// TODO ignore_material datatype
	inline PHYSICS intersect_event propagate(vec3 start, vec3 direction, real distance,
		triangle const * ignore_triangle, int ignore_material) const;

	// Get the maximum distance that can be travelled inside the simulation domain.
	inline PHYSICS real get_max_extent() const;

	// Get the (axis-aligned) simulation domain
	inline PHYSICS vec3 AABB_min() const;
	inline PHYSICS vec3 AABB_max() const;

private:
	CPU void set_AABB(vec3 min, vec3 max);

	triangle* _triangles = nullptr;
	triangle_index_t _N  = 0;
	vec3 _AABB_min       = { 0, 0, 0 };
	vec3 _AABB_max       = { 0, 0, 0 };
	real _max_extent     = 0;

	friend struct _trilist_factory<gpu_flag>;
};

}} // namespace nbl::geometry

#include "trilist.inl"

#endif // __GEOMETRY_TRILIST_H_

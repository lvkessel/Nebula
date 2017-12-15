#include "../legacy_thomas/octree.hh"
#include <stack>

namespace nbl { namespace geometry {

template<bool gpu_flag>
HOST octree<gpu_flag> octree<gpu_flag>::create(std::vector<triangle> const & triangles)
{
	// Conversion between legacy_thomas::point3 and vec3.
	auto p3 = [](vec3 const & v) -> legacy_thomas::point3 { return legacy_thomas::point3(v.x, v.y, v.z); };
	auto v3 = [](legacy_thomas::point3 const & p) -> vec3 { return {(real)p.x, (real)p.y, (real)p.z}; };

	// TODO: error message
	if (triangles.empty())
		throw std::runtime_error("No triangles provided!");

	// Find AABB_min and max
	vec3 AABB_min = triangles[0].AABB_min();
	vec3 AABB_max = triangles[0].AABB_max();

	for (const triangle t : triangles)
	{
		const vec3 tri_min = t.AABB_min();
		const vec3 tri_max = t.AABB_max();

		AABB_min =
		{
			std::min(AABB_min.x, tri_min.x),
			std::min(AABB_min.y, tri_min.y),
			std::min(AABB_min.z, tri_min.z)
		};
		AABB_max =
		{
			std::max(AABB_max.x, tri_max.x),
			std::max(AABB_max.y, tri_max.y),
			std::max(AABB_max.z, tri_max.z)
		};
	}

	AABB_min -= vec3{ 1, 1, 1 };
	AABB_max += vec3{ 1, 1, 1 };

	// Create legacy_thomas::octree structure
	legacy_thomas::octree root(p3(AABB_min), p3(AABB_max));
    for(auto cit = triangles.cbegin(); cit != triangles.cend(); cit++)
	{
		legacy_thomas::triangle legacy_tri(p3(cit->r0()), p3(cit->r1()), p3(cit->r2()), cit->material_in, cit->material_out);
        root.insert(legacy_tri);
	}


	/*
	 * Following is mostly copied from legacy -> cuda_geometry_struct::create()
	 */
	// sort octree nodes by location code (morton order)
    std::map<uint64_t,const legacy_thomas::octree*> morton_map;
    std::stack<const legacy_thomas::octree*> node_p_stack;
    node_p_stack.push(&root);
    while(!node_p_stack.empty()) {
        const legacy_thomas::octree* node_p = node_p_stack.top();
        node_p_stack.pop();
        morton_map[node_p->location()] = node_p;
        for(int octant = 0; octant < 8; octant++) {
            const legacy_thomas::octree* child_p = node_p->traverse(octant);
            if(child_p != nullptr)
                node_p_stack.push(child_p);
        }
    }

    // map triangles from octree to indices following location code
    std::vector<const legacy_thomas::triangle*> triangle_p_vec;
    std::map<const legacy_thomas::triangle*,int> triangle_p_map;
    for(auto morton_cit = morton_map.cbegin(); morton_cit != morton_map.cend(); morton_cit++) {
        const legacy_thomas::octree* node_p = morton_cit->second;
        if(node_p->is_leaf())
            for(const legacy_thomas::triangle* triangle_p : node_p->triangles())
                if(triangle_p_map.count(triangle_p) == 0) {
                    const int index = triangle_p_vec.size();
                    triangle_p_map[triangle_p] = index;
                    triangle_p_vec.push_back(triangle_p);
                }
    }

    // build linearized octree index table
    //  index=0 : child does not exist
    //  index>0 : non-leaf child with node indices
    //  index<0 : leaf child with triangle indices (index -1 means no triangle)
    std::vector<int> octree_vec;
    std::map<const legacy_thomas::octree*,int> node_p_map;
    for(auto morton_cit = morton_map.cbegin(); morton_cit != morton_map.cend(); morton_cit++) {
        const legacy_thomas::octree* node_p = morton_cit->second;
        const int index = octree_vec.size();
        node_p_map[node_p] = index;
        if(node_p->is_leaf()) {
            for(const legacy_thomas::triangle* triangle_p : node_p->triangles())
                octree_vec.push_back(triangle_p_map[triangle_p]);
            octree_vec.push_back(-1);
        } else {
            for(int octant = 0; octant < 8; octant++)
                octree_vec.push_back(0);
        }
    }
    for(auto cit = node_p_map.cbegin(); cit != node_p_map.cend(); cit++) {
        const legacy_thomas::octree* node_p = cit->first;
        const int index = cit->second;
        if(!node_p->is_leaf())
            for(int octant = 0; octant < 8; octant++) {
                const legacy_thomas::octree* child_p = node_p->traverse(octant);
                if(child_p != nullptr) {
                    octree_vec[index+octant] = node_p_map[child_p];
                    if(child_p->is_leaf())
                        octree_vec[index+octant] *= -1;
                }
            }
    }
	/*
	 * END copy from legacy
	 */

	// TODO: ensure that triangles.size() fits in octree._N
	return _octree_factory<gpu_flag>::create(triangles.size(), octree_vec, triangle_p_vec, AABB_min, AABB_max);
}

template<bool gpu_flag>
HOST void octree<gpu_flag>::destroy(octree<gpu_flag> & geometry)
{
	_octree_factory<gpu_flag>::free(geometry);
}

template<bool gpu_flag>
PHYSICS bool octree<gpu_flag>::in_domain(vec3 pos) const
{
	return ((pos.x > _AABB_center.x-_AABB_halfsize.x) && (pos.x < _AABB_center.x+_AABB_halfsize.x) 
		&& (pos.y > _AABB_center.y-_AABB_halfsize.y) && (pos.y < _AABB_center.y+_AABB_halfsize.y)
		&& (pos.z > _AABB_center.z-_AABB_halfsize.z) && (pos.z < _AABB_center.z+_AABB_halfsize.z));
}

template<bool gpu_flag>
PHYSICS intersect_event octree<gpu_flag>::propagate(vec3 start, vec3 direction, real distance,
	triangle const * ignore_triangle, int ignore_material)
{
	// Current (x, y, z) position in the octree, is updated as we travel through cells
	vec3 current_position = start;
	// Distance we still need to travel from current_position, in units of direction
	real distance_to_go = distance;

	// Search in octree
	// A cell in the octree is uniquely identified by a location code.
	// The root cell has location code 1.
	// The location code of a decedant is obtained by shifting the location code
	// to the left by three bits and adding the octant of the decadant to the location code.
	uint64_t location = 1;
	do {
		// define the root axis aligned bounding box
		vec4 AABB{
			_AABB_center.x,
			_AABB_center.y,
			_AABB_center.z,
			1 // scale factor
		};

		// initial index to linearized octree
		int index = 0;

		// traverse to location: we compute the location of the current node and
		//   set the index to the node in the linearized octree array
		// `__clzll` computes number of leading zeros, effectively this is 64-floor(log2(n))
		for (int i = 60 - clz(location); i >= 0; i -= 3) {
			const unsigned int octant = (location >> i) & 7;
			index = _octree_data[index + octant];
			AABB.x += AABB.w*_AABB_halfsize.x*((octant & 1) - 0.5_r);
			AABB.y += AABB.w*_AABB_halfsize.y*(((octant & 2) >> 1) - 0.5_r);
			AABB.z += AABB.w*_AABB_halfsize.z*(((octant & 4) >> 2) - 0.5_r);
			AABB.w *= 0.5_r;
		}

		// traverse to leaf node (which has a linearized octree index smaller than zero)
		while (index >= 0) {
			unsigned int octant = 0;
			octant ^= (current_position.x > AABB.x) ? 1 : 0;
			octant ^= (current_position.y > AABB.y) ? 2 : 0;
			octant ^= (current_position.z > AABB.z) ? 4 : 0;
			location = (location << 3) | octant;
			index = _octree_data[index + octant];
			AABB.x += AABB.w*_AABB_halfsize.x*((octant & 1) - 0.5_r);
			AABB.y += AABB.w*_AABB_halfsize.y*(((octant & 2) >> 1) - 0.5_r);
			AABB.z += AABB.w*_AABB_halfsize.z*(((octant & 4) >> 2) - 0.5_r);
			AABB.w *= 0.5_r;
		}
		// note that a leaf has a guaranteed negative index, the actual index to the leaf
		// node is then obtained by negating.
		index = -index;


		// determine cell/triangle intersections for the leaf
		real intersect; // Distance to next triangle intersection
		int target; // target >= 0 : triangle index
					// target = -1 : cell x intersection
					// target = -2 : cell y intersection
					// target = -4 : cell z intersection

		// find intersection times with each orthogonal plane of the leaf's bounding box,
		// then see which plane is reached first.
		const vec3 t = AABB_intersect(current_position, direction, { AABB.x, AABB.y, AABB.z }, AABB.w*_AABB_halfsize);
		if ((t.x < t.y) && (t.x < t.z)) {
			intersect = t.x;
			target = -1;
		}
		else if ((t.y < t.x) && (t.y < t.z)) {
			intersect = t.y;
			target = -2;
		}
		else {
			intersect = t.z;
			target = -4;
		}


		// iterate triangles in leaf
		while (true) {
			const int triangle_idx = _octree_data[index++];
			if (triangle_idx < 0)
				break;

			// don't intersect with ignore_triangle
			if (_triangles + triangle_idx == ignore_triangle)
				continue;

			// retrieve triangle from global memory
			const triangle this_triangle = _triangles[triangle_idx];

			int mat_idx_out;
			if (dot_product(this_triangle.get_normal(), direction) < 0)
				mat_idx_out = this_triangle.material_in;
			else
				mat_idx_out = this_triangle.material_out;

			// if the outgoing material is the same as current, nothing happens
			// if the triangle represents a detector which can't see the current
			// particle, nothing happens
			if ((mat_idx_out == ignore_material))// || (mat_idx_out == triangle::NOP))
				continue;

			// compute the intersection with the triangle; keep it if it
			// is closer than the current distance to next cell or triangle
			const real t = this_triangle.intersect_ray(current_position, direction);
			if ((t > 0) && (t <= intersect + EPSILON)) {
				intersect = t;
				target = triangle_idx;
			}
		}


		// manage calculated intersections
		if (intersect >= distance_to_go) {
			// EXIT: distance traveled without triangle intersections
			return { distance, nullptr };
		}
		else if (target >= 0) {
			// EXIT: triangle intersection
			return { distance - distance_to_go + intersect, _triangles + target };
		}

		// go to edge of the leaf
		distance_to_go -= intersect;
		current_position += direction * intersect;

		// find adjacent node by bit-magic, and loop
		// the adjacent node is determined solely by location code!
		unsigned int mask = -target;
		unsigned int value;
		// mask = 1 (001), value = 0 (000) : xx0 --> xx1 (find neighbor in positive x direction)
		// mask = 1 (001), value = 1 (000) : xx1 --> xx0 (find neighbor in negative x direction)
		// mask = 2 (010), value = 0 (000) : x0x --> x1x (find neighbor in positive y direction)
		// mask = 2 (010), value = 2 (010) : x1x --> x0x (find neighbor in negative y direction)
		// mask = 4 (100), value = 0 (000) : 0xx --> 1xx (find neighbor in positive z direction)
		// mask = 4 (100), value = 4 (100) : 1xx --> 0xx (find neighbor in negative z direction)
		if (mask == 1)
			value = (direction.x >= 0) ? 0 : 1;
		else if (mask == 2)
			value = (direction.y >= 0) ? 0 : 2;
		else
			value = (direction.z >= 0) ? 0 : 4;
		while (location > 1) {
			if ((location&mask) == value) {
				location ^= mask;
				break;
			}
			location >>= 3;
		}

	} while (location > 1);

	// EXIT: particle is out of grid (left root cell)
	return { distance, nullptr };
}

template<bool gpu_flag>
PHYSICS real octree<gpu_flag>::get_max_extent() const
{
	return _max_extent;
}

template<bool gpu_flag>
inline PHYSICS vec3 octree<gpu_flag>::AABB_min() const
{
	return _AABB_center - _AABB_halfsize;
}
template<bool gpu_flag>
inline PHYSICS vec3 octree<gpu_flag>::AABB_max() const
{
	return _AABB_center + _AABB_halfsize;
}

template<bool gpu_flag>
HOST void octree<gpu_flag>::set_AABB(vec3 min, vec3 max)
{
	_AABB_center = (min+max)/2;
	_AABB_halfsize.x = std::fabs(max.x-min.x)/2;
	_AABB_halfsize.y = std::fabs(max.y-min.y)/2;
	_AABB_halfsize.z = std::fabs(max.z-min.z)/2;

	const vec3 m = max - min;
	_max_extent = magnitude(m);
}

template<bool gpu_flag>
PHYSICS vec3 octree<gpu_flag>::AABB_intersect(vec3 pos, vec3 dir, vec3 center, vec3 halfsize)
{
	return
	{
		(center.x + copysignr(halfsize.x + EPSILON, dir.x) - pos.x) / dir.x,
		(center.y + copysignr(halfsize.y + EPSILON, dir.y) - pos.y) / dir.y,
		(center.z + copysignr(halfsize.z + EPSILON, dir.z) - pos.z) / dir.z
	};
}

#ifdef _MSC_VER
#include <intrin.h>
#endif

template<bool gpu_flag>
PHYSICS int octree<gpu_flag>::clz(uint64_t x)
{
	static_assert(sizeof(x) == sizeof(long long int), "__clzll intrinsic should be 64 bit?");
#if defined __CUDA_ARCH__
	return __clzll(x);
#elif defined __GNUC__
	return __builtin_clzll(x);
#elif defined(_MSC_VER)
	return __lzcnt64(x);
#else
	static_assert(0, "No clz implementation for your platform!");
#endif
}

template<>
struct _octree_factory<false>
{
	inline static HOST octree<false> create(octree<false>::triangle_index_t N, std::vector<int> octree_vec,
		std::vector<const legacy_thomas::triangle*> triangle_p_vec, vec3 AABB_min, vec3 AABB_max)
	{
		auto v3 = [](legacy_thomas::point3 const & p) -> vec3 { return { (real)p.x, (real)p.y, (real)p.z }; };

		octree<false> geometry;

		geometry._N = N;

		geometry._octree_data = new int[octree_vec.size()];
		memcpy(geometry._octree_data, octree_vec.data(), octree_vec.size() * sizeof(int));

		geometry._triangles = reinterpret_cast<triangle*>(malloc(geometry._N * sizeof(triangle)));
		for (octree<false>::triangle_index_t i = 0; i < triangle_p_vec.size(); ++i)
		{
			legacy_thomas::triangle const * legacy_tri = triangle_p_vec[i];
			geometry._triangles[i] = triangle(v3(legacy_tri->A), v3(legacy_tri->B), v3(legacy_tri->C), legacy_tri->in, legacy_tri->out);
		}

		geometry.set_AABB(AABB_min, AABB_max);

		return geometry;
	}

	inline static HOST void free(octree<false> & geometry)
	{
		delete[] geometry._octree_data;
		::free(geometry._triangles);

		geometry._octree_data = nullptr;
		geometry._triangles = nullptr;
		geometry._N = 0;
	}
};

#if CUDA_AVAILABLE
template<>
struct _octree_factory<true>
{
	inline static HOST octree<true> create(octree<true>::triangle_index_t N, std::vector<int> octree_vec,
		std::vector<const legacy_thomas::triangle*> triangle_p_vec, vec3 AABB_min, vec3 AABB_max)
	{
		auto v3 = [](legacy_thomas::point3 const & p) -> vec3 { return { (real)p.x, (real)p.y, (real)p.z }; };

		octree<true> geometry;

		geometry._N = N;

		// Copy octree to device
		cuda::cuda_new<int>(&geometry._octree_data, octree_vec.size());
		cudaMemcpy(geometry._octree_data, octree_vec.data(), octree_vec.size() * sizeof(int), cudaMemcpyHostToDevice);

		// Copy triangle data to device
		cuda::cuda_new<triangle>(&geometry._triangles, geometry._N);
		cuda::cuda_mem_scope<triangle>(geometry._triangles, geometry._N,
			[&triangle_p_vec, &v3](triangle* device)
		{
			for (octree<true>::triangle_index_t i = 0; i < triangle_p_vec.size(); ++i)
			{
				legacy_thomas::triangle const * legacy_tri = triangle_p_vec[i];
				device[i] = triangle(v3(legacy_tri->A), v3(legacy_tri->B), v3(legacy_tri->C), legacy_tri->in, legacy_tri->out);
			}
		});

		geometry.set_AABB(AABB_min, AABB_max);

		return geometry;
	}

	inline static HOST void free(octree<true> & geometry)
	{
		cudaFree(geometry._octree_data);
		cudaFree(geometry._triangles);

		geometry._octree_data = nullptr;
		geometry._triangles = nullptr;
		geometry._N = 0;
	}
};
#endif // CUDA_AVAILABLE

}} // namespace nbl::geometry

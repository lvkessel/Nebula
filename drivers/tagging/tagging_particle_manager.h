#ifndef __TAGGING_PARTICLE_MANAGER_H_
#define __TAGGING_PARTICLE_MANAGER_H_

#include <vector>

namespace nbl { namespace drivers {

template<typename material_manager_t, typename cascade_info>
class tagging_particle_manager
{
public:
	using status_t = uint8_t;
	using particle_index_t = size_t;
	using material_index_t = typename material_manager_t::material_index_t;
	using primary_tag_t = uint32_t;
	using unique_tag_t = uint64_t;

	static tagging_particle_manager create();
	static void destroy(tagging_particle_manager & manager);

	// Push particles to GPU. Returns how many particles were actually pushed.
	particle_index_t push(particle* particles, primary_tag_t* tags, particle_index_t N);

	template<typename detect_function>
	void flush_detected(detect_function func);
	void flush_terminated();

	particle_index_t get_total_count() const;
	particle_index_t get_running_count() const;
	particle_index_t get_detected_count() const;

	// Direct access, no bounds checking
	inline PHYSICS particle & operator[](particle_index_t i);
	inline PHYSICS particle const & operator[](particle_index_t i) const;

	// Is there a memory location for this particle?
	inline PHYSICS bool exists(particle_index_t i) const;
	// Is particle active == not PENDING, DETECTED or TERMINATED
	inline PHYSICS bool active(particle_index_t i) const;

	// Get current material
	inline PHYSICS material_index_t get_material_index(particle_index_t i) const;
	inline PHYSICS void set_material_index(particle_index_t particle_idx, material_index_t new_material_idx);

	// Get last intersected triangle for a particle (or nullptr)
	inline PHYSICS triangle const * get_last_triangle(particle_index_t i) const;
	inline PHYSICS void forget_last_triangle(particle_index_t i);

	// Is next event scatter / intersect?
	inline PHYSICS bool next_scatter(particle_index_t i) const;
	inline PHYSICS uint8_t get_next_scatter(particle_index_t i) const;
	inline PHYSICS bool next_intersect(particle_index_t i) const;

	// Add new particle
	inline PHYSICS void create_secondary(particle_index_t primary_idx, particle secondary_particle);

	// Terminate particle
	inline PHYSICS void terminate(particle_index_t i);
	// Detect particle
	inline PHYSICS void detect(particle_index_t i);
	// Set next scattering event
	inline PHYSICS void set_scatter_event(particle_index_t i, scatter_event event);
	inline PHYSICS void set_intersect_event(particle_index_t i, intersect_event event);

protected:
	struct particle_struct : cascade_info
	{
		particle_struct(status_t _status, uint8_t _next_scatter, material_index_t _current_material, particle _particle_data,
			primary_tag_t _primary_tag, unique_tag_t _unique_tag, triangle* _last_triangle, cascade_info _info = cascade_info{})
			: cascade_info(_info), status(_status), next_scatter(_next_scatter), current_material(_current_material),
			particle_data(_particle_data), primary_tag(_primary_tag), unique_tag(_unique_tag), last_triangle(_last_triangle)
		{}

		status_t status;
		uint8_t next_scatter;
		material_index_t current_material;
		particle particle_data;
		primary_tag_t primary_tag;
		unique_tag_t unique_tag;
		triangle* last_triangle;
	};

	std::vector<particle_struct> data;
	unique_tag_t next_unique_tag = 0;

	enum status_enum : status_t
	{
		SCATTER_EVENT   = 0b0100,
		INTERSECT_EVENT = 0b0010,
		NO_EVENT        = 0b0110,
		DETECTED        = 0b1010,
		NEW_SECONDARY   = 0b1110,
		TERMINATED      = 0b0011
	};
};

class empty_class {};
template<typename material_manager_t>
using default_tagging_particle_manager = tagging_particle_manager<material_manager_t, empty_class>;

}} // namespace nbl::drivers

#include "tagging_particle_manager.inl"

#endif // __TAGGING_PARTICLE_MANAGER_H_

#ifndef __CPU_PARTICLE_MANAGER_H_
#define __CPU_PARTICLE_MANAGER_H_

#include <vector>

namespace nbl { namespace drivers {

/*
 * Generic CPU particle manager.
 * Serves as a base class for a "simple" CPU particle manager
 * and a version with extended tracking facilities.
 */
template<typename material_manager_t, typename additional_data>
class cpu_particle_manager
{
public:
	using status_t = uint8_t;
	using particle_index_t = size_t;
	using material_index_t = typename material_manager_t::material_index_t;
	using primary_tag_t = uint32_t;

	static cpu_particle_manager create();
	static void destroy(cpu_particle_manager & manager);

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
	/*
	 * particle_struct holds relevant data for each particle.
	 * It has a template parameter, T, which is inherited from if not void.
	 * This can be used to store additional data.
	 *
	 * The int template parameter is there to make sure that particle_struct<void>
	 * is a PARTIAL specialization. For some reason, that's allowed, while a full
	 * specialization is not.
	 */
	template<typename T, int = 0>
	struct particle_struct : T
	{
		status_t status;
		uint8_t next_scatter;
		material_index_t current_material;
		particle particle_data;
		primary_tag_t primary_tag;
		triangle* last_triangle;

		// Provide constructor.
		// P0017R1 (C++17) makes this superfluous.
		particle_struct(T const & extended_data,
			status_t status, uint8_t next_scatter, material_index_t current_material,
			particle particle_data, primary_tag_t primary_tag, triangle* last_triangle)
			: T(extended_data), status(status), next_scatter(next_scatter),
			current_material(current_material), particle_data(particle_data),
			primary_tag(primary_tag), last_triangle(last_triangle)
		{}
	};
	template<int i>
	struct particle_struct<void, i>
	{
		status_t status;
		uint8_t next_scatter;
		material_index_t current_material;
		particle particle_data;
		primary_tag_t primary_tag;
		triangle* last_triangle;
	};

	std::vector<particle_struct<additional_data>> data;

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

}} // namespace nbl::drivers

#include "cpu_particle_manager.inl"

#endif // __CPU_PARTICLE_MANAGER_H_

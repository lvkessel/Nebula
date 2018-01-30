#ifndef __CPU_DRIVER_H_
#define __CPU_DRIVER_H_

#include "../../core/material_manager.h"
#include "../../common/util/random.h"
#include "simple_cpu_particle_manager.h"

namespace nbl { namespace drivers {

template<
	typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager = simple_cpu_particle_manager
>
class cpu_driver
{
public:
	using material_t = material<scatter_list_t>;
	using material_manager_t = material_manager<material_t, false>;
	using particle_manager_t = particle_manager<material_manager_t>;

	using particle_index_t = typename particle_manager_t::particle_index_t;
	using primary_tag_t = typename particle_manager_t::primary_tag_t;
	using status_t = typename particle_manager_t::status_t;

	cpu_driver(geometry_manager_t geometry, intersect_t intersect, std::vector<material_t> materials);
	~cpu_driver();

	inline particle_index_t push(particle* particles, primary_tag_t* tags, particle_index_t N);
	void do_iteration();
	void simulate_until_end();

	particle_index_t get_running_count() const;
	particle_index_t get_detected_count() const;

	// Detect_function is called for each particle before it is terminated.
	// One will typically pass a function that writes interesting data to an output stream.
	template<typename detect_function>
	void flush_detected(detect_function function);

	void flush_terminated();

private:
	particle_manager_t _particles;
	material_manager_t _materials;
	geometry_manager_t _geometry;
	intersect_t _intersect;

	// Random number generators
	util::random_generator<false> rand_state;

	// Functions called by do_iteration()
	inline void init(particle_index_t particle_idx);
	inline void intersect(particle_index_t particle_idx);
	inline void scatter(particle_index_t particle_idx);

	cpu_driver(cpu_driver const &) = delete;
	cpu_driver& operator=(cpu_driver const &) = delete;
};

}} // namespace nbl::drivers

#include "cpu_driver.inl"

#endif // __CPU_DRIVER_H_

#ifndef __GPU_DRIVER_H_
#define __GPU_DRIVER_H_

/*
 * Driver similar to Thomas Verduin's e-scatter.
 */

#include "../../core/material_manager.h"
#include "../../common/util/random.h"
#include "gpu_particle_manager.h"

namespace nbl { namespace drivers {

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
class gpu_driver
{
public:
	static_assert(scatter_list_t::size() == 2,
		"GPU driver must have two scattering mechanisms");
	static_assert(scatter_list_t::template type_at_index<1>::may_create_se == false,
		"Only the first scattering mechanism may create secondaries in GPU code");

	using material_t = material<scatter_list_t>;
	using material_manager_t = material_manager<material_t, true>;
	using particle_manager_t = gpu_particle_manager<material_manager_t>;

	using particle_index_t = typename particle_manager_t::particle_index_t;
	using status_t = typename particle_manager_t::status_t;

	HOST gpu_driver(particle_index_t particle_capacity, geometry_manager_t geom,
		intersect_t inter, std::vector<material_t> materials);
	HOST ~gpu_driver();

	inline HOST particle_index_t push(particle* particles, uint32_t* tags, particle_index_t N);
	HOST void do_iteration();

	HOST particle_index_t get_running_count() const;
	HOST particle_index_t get_detected_count() const;
	
	// Detect_function is called for each particle before it is terminated.
	// One will typically pass a function that writes interesting data to an output stream.
	template<typename detect_function>
	HOST void flush_detected(detect_function function);
	
private:
	particle_manager_t _particles;
	material_manager_t _materials;
	geometry_manager_t _geometry;
	intersect_t _intersect;

	// GPU run settings
	unsigned int _threads_per_block = 32;
	unsigned int _num_blocks = 0;

	// Random number generators
	util::random_generator<true>* curand_states = nullptr;

	// Functions called by do_iteration()
	inline HOST void init();
	inline HOST void events();

	gpu_driver(gpu_driver const &) = delete;
	gpu_driver& operator=(gpu_driver const &) = delete;
};

}} // namespace nbl::drivers

#include "gpu_driver.inl"

#endif // __GPU_DRIVER_H_

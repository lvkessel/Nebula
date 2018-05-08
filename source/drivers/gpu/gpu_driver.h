#ifndef __GPU_DRIVER_H_
#define __GPU_DRIVER_H_

#include "../../core/material_manager.h"
#include "../../common/util/random.h"
#include "gpu_particle_manager.h"

namespace nbl { namespace drivers {

/**
 * \brief Driver similar to Thomas Verduin's e-scatter.
 *
 * It only accepts two scattering mechanisms. The first of these may generate
 * secondary electrons, the second may not. Only one SE may be generated per
 * event.
 *
 * T.V.'s thesis: doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c
 *
 * \tparam scatter_list_t     List of scattering mechanisms. Should be of type
 *                            ::scatter_list.
 * \tparam intersect_t        Intersection event handler
 * \tparam geometry_manager_t Geometry manager
 */
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
	using seed_t = typename util::random_generator<true>::seed_type;

	/**
	 * \brief Constructor.
	 *
	 * \param particle_capacity Maximum number of particles the simulation may have.
	 * \param geom              Geometry manager, holding the simulation geometry.
	 * \param inter             Instance of the intersection handler.
	 * \param materials         List of materials in the simulation.
	 * \param seed              Seed for the random number generator.
	 */
	CPU gpu_driver(particle_index_t particle_capacity, geometry_manager_t geom,
		intersect_t inter, std::vector<material_t> materials,
		seed_t seed = util::random_generator<true>::default_seed);
	CPU ~gpu_driver();

	/**
	 * \brief Add new primary particles to the simulation.
	 *
	 * \param particles Pointer to an array of particles to be added.
	 * \param tags      Pointer to the corresponding array of tags.
	 * \param N         Number of particles to be added.
	 */
	inline CPU particle_index_t push(particle* particles, uint32_t* tags, particle_index_t N);

	/// Perform a single iteration of the simulation for all particles.
	CPU void do_iteration();

	/**
	 * \brief Get number of particles currently in the simulation.
	 *
	 * More specifically, those that are not terminated or detected.
	 */
	CPU particle_index_t get_running_count() const;

	/// Get number of detected particles currently in the simulation.
	CPU particle_index_t get_detected_count() const;

	/**
	 * \brief Set the detected particles to terminated, calling a callback
	 *        before doing so.
	 *
	 * The callback function receives the particle data and the tag that
	 * belonged to the primary electron that initated the cascade the detected
	 * particle is part of.
	 *
	 * \param function Callback function to be called for each detected particle.
	 *                 Should have signature void(particle const &, uint32_t).
	 */
	template<typename detect_function>
	CPU void flush_detected(detect_function function);

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
	inline CPU void init();
	inline CPU void events();

	gpu_driver(gpu_driver const &) = delete;
	gpu_driver& operator=(gpu_driver const &) = delete;
};

}} // namespace nbl::drivers

#include "gpu_driver.inl"

#endif // __GPU_DRIVER_H_

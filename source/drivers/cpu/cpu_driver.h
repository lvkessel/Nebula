#ifndef __CPU_DRIVER_H_
#define __CPU_DRIVER_H_

#include "../../core/material_manager.h"
#include "../../common/util/random.h"
#include "simple_cpu_particle_manager.h"

namespace nbl { namespace drivers {

/**
 * \brief Generic CPU driver
 *
 * Accepts any combination of scattering mechanisms, any number of secondaries
 * may be generated per event.
 *
 * The particle manager can be provided as a template parameter. The default
 * choice is perfect for most simulations, but users may want to use custom
 * particle managers if, for example, they want to track the cascade in more
 * detail.
 *
 * \tparam scatter_list_t     List of scattering mechanisms. Should be of type
 *                            ::scatter_list.
 * \tparam intersect_t        Intersection event handler
 * \tparam geometry_manager_t Geometry manager
 * \tparam particle_manager   Particle manager to be used.
 */
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
	using seed_t = typename util::random_generator<false>::seed_type;

	/**
	 * \brief Constructor.
	 *
	 * \param geometry  Geometry manager, holding the simulation geometry.
	 * \param intersect Instance of the intersection handler.
	 * \param materials List of materials in the simulation.
	 * \param seed      Seed for the random number generator.
	 */
	cpu_driver(geometry_manager_t geometry, intersect_t intersect, std::vector<material_t> materials,
		seed_t seed = util::random_generator<false>::default_seed);
	~cpu_driver();

	/**
	 * \brief Add new primary particles to the simulation.
	 *
	 * \param particles Pointer to an array of particles to be added.
	 * \param tags      Pointer to the corresponding array of tags.
	 * \param N         Number of particles to be added.
	 */
	inline particle_index_t push(particle* particles, primary_tag_t* tags, particle_index_t N);

	/// Perform a single iteration of the simulation for all particles.
	void do_iteration();

	/// Keep simulating until there are no particles left.
	void simulate_until_end();

	/**
	 * \brief Get number of particles currently in the simulation.
	 *
	 * More specifically, those that are not terminated or detected.
	 */
	particle_index_t get_running_count() const;

	/// Get number of detected particles currently in the simulation.
	particle_index_t get_detected_count() const;

	/**
	 * \brief Set the detected particles to terminated, calling a callback
	 *        before doing so.
	 *
	 * The callback function receives the particle data and the tag that
	 * belonged to the primary electron that initated the cascade the detected
	 * particle is part of.
	 *
	 * These detected particles are simply set to "terminated", meaning they are
	 * still in memory.
	 *
	 * \param function Callback function to be called for each detected particle.
	 *                 Should have signature void(particle const &, uint32_t).
	 */
	template<typename detect_function>
	void flush_detected(detect_function function);

	/// Remove all terminated particles from memory.
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

#ifndef __WORK_POOL_H_
#define __WORK_POOL_H_

#include <mutex>
#include "../core/particle.h"

/**
 * \brief Very simple class to keep track of work to be done in a thread-safe way.
 *
 * This class does not take ownership of the data.
 * It simply assumes that the particles are not deleted until the work is done.
 */

class work_pool
{
public:
	/**
	 * \brief Constructor.
	 *
	 * \param primaries Pointer to primary electron data
	 * \param tags      Pointer to tag data
	 * \param N         Number of particles to be simulated
	 */
	work_pool(particle* primaries, uint32_t* tags, size_t N) :
		next_primary(primaries), next_tag(tags), primaries_to_go(N)
	{}

	/**
	 * \brief Get work to be done
	 *
	 * \param batch_size Number of particles requested
	 *
	 * \return Tuple of
	 *   - Pointer to first particle
	 *   - Pointer to first corresponding tag
	 *   - Number of particles obtained, may be less than batch_size
	 */
	std::tuple<particle*, uint32_t*, size_t> get_work(size_t batch_size)
	{
		std::lock_guard<std::mutex> lock(mutex);

		auto particles_pushed = std::min(batch_size, primaries_to_go);
		auto return_data = std::make_tuple(next_primary, next_tag, particles_pushed);
		next_primary += particles_pushed;
		next_tag += particles_pushed;
		primaries_to_go -= particles_pushed;

		return return_data;
	}

	/// Get the amount of work left in the pool
	size_t get_primaries_to_go() const
	{
		std::lock_guard<std::mutex> lock(mutex);
		return primaries_to_go;
	}

	/// Test if the amount of work is equal to zero.
	bool done() const
	{
		return get_primaries_to_go() == 0;
	}

private:
	mutable std::mutex mutex;

	particle* next_primary;
	uint32_t* next_tag;
	size_t primaries_to_go;
};

#endif // __WORK_POOL_H_

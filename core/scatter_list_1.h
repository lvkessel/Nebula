#ifndef __SCATTER_LIST_H_
#define __SCATTER_LIST_H_

/*
 * scatter_list represents a list of scattering processes
 * that can take place in a material.
 * 
 * It exposes two functions:
 *   sample_path() returns a scatter_event describing which event
 *                 takes place first, and at what distance
 *   execute() executes the desired event.
 * 
 * Each scatter_type is a class that must have public member functions
 *   PHYSICS real sample_path() const;
 *   PHYSICS void execute_event(particle&) const;
 */

#include "tuple.h"
#include "events.h"
#include "particle.h"

template<typename... scatter_types>
class scatter_list : nbl::tuple::tuple<scatter_types...>
{
public:
	HOST scatter_list(scatter_types... sc)
		: nbl::tuple::tuple<scatter_types...>(sc...)
	{}

	// TODO: rng_t --> gpu_flag
	template<typename rng_t>
	inline PHYSICS scatter_event sample_path(particle const & this_particle, rng_t & rng) const
	{
		scatter_event evt{ 0, std::numeric_limits<real>::infinity() };
		nbl::tuple::visit(*this, [&](auto scatter, size_t i)
			{
				auto path = scatter.sample_path(this_particle, rng);
				if (path < evt.distance)
					evt = { uint8_t(2-i), path };
			});

		return evt;
	}

	// TODO: rng_t --> gpu_flag
	template<typename particle_manager, typename rng_t>
	inline PHYSICS void execute(uint8_t i,
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		rng_t& rng) const
	{
		nbl::tuple::visit_at(*this, 2-i,
			[&](auto scatter)
			{ scatter.execute(particle_mgr, particle_idx, rng); });
	}
};

#endif // __SCATTER_LIST_H_
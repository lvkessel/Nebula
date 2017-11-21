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

// TODO namespace

template<typename... scatter_types>
class scatter_list : nbl::tuple::tuple<scatter_types...>
{
public:
	HOST scatter_list(scatter_types... sc)
		: nbl::tuple::tuple<scatter_types...>(sc...)
	{}

	template <size_t I>
	using type_at_index = typename nbl::tuple::type_at_index<I, scatter_types...>::type;

	static constexpr size_t size() { return sizeof...(scatter_types); }

/*
	// Generic lambdas... would be nice, right?
	// But we chose not to rely on C++14 yet...

	template<typename rng_t>
	inline PHYSICS scatter_event sample_path(particle const & this_particle, rng_t & rng) const
	{

		scatter_event evt{ 0, std::numeric_limits<real>::infinity() };
		nbl::tuple::visit(*this, [&](auto scatter, size_t i)
			{
				auto path = scatter.sample_path(this_particle, rng);
				if (path < evt.distance)
					evt = { uint8_t(i+1), path };
			});
		return evt;
	}

	template<typename particle_manager, typename rng_t>
	inline PHYSICS void execute(uint8_t i,
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		rng_t& rng) const
	{
		nbl::tuple::visit_at(*this, i-1,
			[&](auto scatter)
			{ scatter.execute(particle_mgr, particle_idx, rng); });
	}
*/

	// C++11 workarounds
	template<typename rng_t>
	inline PHYSICS scatter_event sample_path(particle const & this_particle, rng_t & rng) const
	{
		sample_visitor<rng_t> visitor{ this_particle, rng, { 0, std::numeric_limits<real>::infinity() }};
		nbl::tuple::visit(*this, visitor);
		return visitor.evt;
	}
	template<typename particle_manager, typename rng_t>
	inline PHYSICS void execute(uint8_t i,
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		rng_t& rng) const
	{
		execute_visitor<particle_manager, rng_t> visitor{ particle_mgr, particle_idx, rng };
		nbl::tuple::visit_at(*this, i-1, visitor);
	}

private:
	template<typename rng_t>
	struct sample_visitor
	{
		particle const & this_particle;
		rng_t & rng;
		scatter_event evt;

		template<typename scatter_type>
		inline PHYSICS void operator()(scatter_type& scatter, size_t i)
		{
			auto path = scatter.sample_path(this_particle, rng);
			if (path < evt.distance)
				evt = { uint8_t(i+1), path };
		}
	};
	template<typename particle_manager, typename rng_t>
	struct execute_visitor
	{
		particle_manager& particle_mgr;
		typename particle_manager::particle_index_t& particle_idx;
		rng_t& rng;

		template<typename scatter_type>
		inline PHYSICS void operator()(scatter_type& scatter)
		{
			scatter.execute(particle_mgr, particle_idx, rng);
		}
	};
};

#endif // __SCATTER_LIST_H_

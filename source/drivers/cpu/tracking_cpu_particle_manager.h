#ifndef __TRACKING_CPU_PARTICLE_MANAGER_H_
#define __TRACKING_CPU_PARTICLE_MANAGER_H_

#include "cpu_particle_manager.h"
#include <tuple>

namespace nbl { namespace drivers {


// This is the data to be stored for every iteration
struct event_info
{
	uint8_t type;  // Event type. Zero is no event, 1 and later correspond to the order in the material manager.
	vec3 position; // Position at which the event took place
	real energy;   // Energy just before scattering
};

// This is the additional data to be stored for every electron
struct cascade_info
{
	uint64_t unique_tag;            // Unique tag belonging to this particle
	uint64_t parent_unique_tag;     // Unique tag that belonged to the parent
	size_t parent_create_event;     // Index to the parent's "events" vector in which this electron was created
	std::vector<event_info> events; // Vector of all scattering events
};

template<typename material_manager_t>
class tracking_cpu_particle_manager
	: public cpu_particle_manager<material_manager_t, cascade_info>
{
public:
	using base_t = cpu_particle_manager<material_manager_t, cascade_info>;
	using typename base_t::particle_index_t;
	using typename base_t::primary_tag_t;
	using base_t::data;

	using track_iterator = typename std::vector<typename base_t::particle_struct>::iterator;

	tracking_cpu_particle_manager(base_t const & base)
		: base_t(base)
	{}

	static tracking_cpu_particle_manager create()
	{
		tracking_cpu_particle_manager<material_manager_t> manager(base_t::create());
		return manager;
	}
	static void destroy(tracking_cpu_particle_manager & manager)
	{
		base_t::destroy(manager);
	}

	particle_index_t push(particle* particles, primary_tag_t* tags, particle_index_t N)
	{
		data.reserve(data.size() + N);
		for (particle_index_t i = 0; i < N; ++i)
		{
			data.push_back({
				base_t::NO_EVENT,
				0,
				-123,  // TODO: vacuum
				particles[i],
				tags[i],
				nullptr,
				{++next_unique_tag, 0, 0, {}}
			});
		}
		return N;
	}

	inline PHYSICS void create_secondary(particle_index_t primary_idx, particle secondary_particle)
	{
		data.push_back({
			base_t::NEW_SECONDARY,
			0,
			base_t::get_material_index(primary_idx),
			secondary_particle,
			data[primary_idx].primary_tag,
			nullptr,
			{++next_unique_tag, data[primary_idx].unique_tag, data[primary_idx].events.size()-1, {}}
		});
	}

	// We want to store all scattering events
	inline PHYSICS void set_scatter_event(particle_index_t i, scatter_event event)
	{
		// Base class's set_scatter_event moves the particle's position
		base_t::set_scatter_event(i, event);
		data[i].events.push_back({ event.type, data[i].particle_data.pos, data[i].particle_data.kin_energy });
	}

	/*
	 * Overrride the standard detection function.
	 * We leave the "detection function" the same -- the main simulator just gets
	 * the particle data and primary tag for each detected electron.
	 * In addition, we will, for each detected electron, write the deepest point,
	 * energy at deepest point and energy at detection to stdout. We consider the
	 * electron itself as well as all parents.
	 */
	template<typename detect_function>
	void flush_detected(detect_function func)
	{
		for (auto it = data.begin(); it != data.end(); ++it)
		{
			if (it->status == base_t::DETECTED)
			{
				// Execute the regular "detection function" and terminate the particle
				func(it->particle_data, it->primary_tag);
				it->status = base_t::TERMINATED;


				// This is what we want to know.
				real deepest_z = std::numeric_limits<real>::infinity();
				real deepest_energy = -1;

				/*
				 * Loop through all events that this particle and its parents have had.
				 * For a parent, we are only interested in all events before final_event_index,
				 * which is the event at which a parent created the child.
				 */
				auto current_parent = it;
				auto final_event_index = it->events.size();
				// Loop through all parents
				while (current_parent != data.end())
				{
					// Loop through all events before the child was created,
					// looking for an event deeper than the one we already have
					for (size_t event_index = 0; event_index < final_event_index; ++event_index)
					{
						const event_info ei = current_parent->events[event_index];
						if (ei.position.z < deepest_z)
						{
							deepest_z = ei.position.z;
							deepest_energy = ei.energy;
						}
					}

					// Go to next parent
					final_event_index = current_parent->parent_create_event;
					current_parent = find_parent(current_parent);
				}

				// Send data to stdout.
				// It is possible that this particle was detected without any events,
				// in which case deepest_energy is still -1. Do not print anything
				// in that case.
				if (deepest_energy > 0)
				{
					std::cout << deepest_z << '\t'
						<< deepest_energy << '\t'
						<< it->particle_data.kin_energy << '\n';
				}
			}
		}
	}


#if 0
	/*
	 * Alternative version, that does not check all parents.
	 * (Disabled by the #if 0 above and #endif below)
	 */
	template<typename detect_function>
	void flush_detected(detect_function func)
	{
		for (auto& this_particle : data)
		{
			if (this_particle.status == base_t::DETECTED)
			{
				// Execute the regular "detection function" and terminate the particle
				func(this_particle.particle_data, this_particle.primary_tag);
				this_particle.status = base_t::TERMINATED;

				// Find the deepest point of this particle (but not its parents)
				real deepest_z = std::numeric_limits<real>::infinity();
				real deepest_energy = -1;
				for (event_info ei : this_particle.events)
				{
					if (ei.position.z < deepest_z)
					{
						deepest_z = ei.position.z;
						deepest_energy = ei.energy;
					}
				}

				// And output interesting data
				if (deepest_energy > 0)
				{
					std::cout << deepest_z << '\t'
						<< deepest_energy << '\t'
						<< this_particle.particle_data.kin_energy << '\n';
				}
			}
		}
	}
#endif

protected:
	/*
	 * Find a particle's track data by its unique tag.
	 * A unique tag stays the same throughout the simulation, while these
	 * iterators may be invalidated.
	 * Returns data.end() if not found.
	 */
	track_iterator find_by_unique_tag(uint64_t unique_tag)
	{
		// Note that, while particle's indices may change during a simulation,
		// their order does not. They are sorted by unique tag.
		// We use this to speed up the search.
		auto iterator = std::lower_bound(data.begin(), data.end(), unique_tag,
			[](cascade_info const & a, uint64_t b) -> bool
			{ return a.unique_tag < b; });

		if (iterator->unique_tag != unique_tag)
			return data.end();

		return iterator;
	}

	/*
	 * Find a particle's parent.
	 * Slightly faster than find_by_unique_tag(iterator->parent_unique_tag),
	 * because we don't need to search the whole range.
	 */
	track_iterator find_parent(track_iterator child)
	{
		auto iterator = std::lower_bound(data.begin(), child, child->parent_unique_tag,
			[](cascade_info const & a, uint64_t b) -> bool
			{ return a.unique_tag < b; });

		if (iterator->unique_tag != child->parent_unique_tag)
			return data.end();

		return iterator;
	}

	uint64_t next_unique_tag = 0;
};

}} // namespace nbl::drivers

#endif // __TRACKING_CPU_PARTICLE_MANAGER_H_

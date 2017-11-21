#ifndef __CASCADE_SAVING_PARTICLE_MANAGER_H_
#define __CASCADE_SAVING_PARTICLE_MANAGER_H_

#include "tagging_particle_manager.h"
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
	// TODO: data type
	uint64_t parent_unique_tag = 0; // Unique tag that belonged to the parent
	size_t parent_create_event = 0; // Index to the parent's "events" vector in which this electron was created
	std::vector<event_info> events; // Vector of all scattering events
};

template<typename material_manager_t>
class cascade_saving_particle_manager
	: public tagging_particle_manager<material_manager_t, cascade_info>
{
public:
	using base_t = tagging_particle_manager<material_manager_t, cascade_info>;
	using typename base_t::particle_index_t;
	using base_t::data;

	cascade_saving_particle_manager(base_t const & base)
		: base_t(base)
	{}

	static cascade_saving_particle_manager create()
	{
		cascade_saving_particle_manager<material_manager_t> manager(base_t::create());
		return manager;
	}
	static void destroy(cascade_saving_particle_manager & manager)
	{
		base_t::destroy(manager);
	}

	// No need to override push: default constructor for cascade_info is OK.
	// Do need to override create_secondary: need to provide appropriate parent tag and event # of creation
	inline PHYSICS void create_secondary(particle_index_t primary_idx, particle secondary_particle)
	{
		base_t::create_secondary(primary_idx, secondary_particle);
		data.back().parent_unique_tag = data[primary_idx].unique_tag;
		data.back().parent_create_event = data[primary_idx].events.size() - 1;
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
	 * In addition, this particle manager will, for each detected electron,
	 * write the deepest point, energy at deepest point and energy at detection to stdout.
	 * It considers the electron itself as well as all its parents.
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


				// This is what we want to know.
				real deepest_z = std::numeric_limits<real>::infinity();
				real deepest_energy = -1;

				/*
				 * Loop through all events that this particle and its parents have had.
				 * The "current parent" is called current_particle here (sorry for bad naming).
				 * For a parent, we are only interested in all events before final_event_index,
				 * which is the event at which a parent created the child.
				 */
				auto const * current_particle = &this_particle;
				auto final_event_index = this_particle.events.size();
				while (current_particle != nullptr)
				{
					// See if this particle has had an event even deeper than what we already have
					for (size_t event_index = 0; event_index < final_event_index; ++event_index)
					{
						const event_info ei = current_particle->events[event_index];
						if (ei.position.z < deepest_z)
						{
							deepest_z = ei.position.z;
							deepest_energy = ei.energy;
						}
					}


					// Go to next parent
					final_event_index = current_particle->parent_create_event;
					const auto parent_unique_tag = current_particle->parent_unique_tag;
					current_particle = nullptr;
					if (parent_unique_tag != 0) // == 0 means primary particle
					{
						for (auto& candidate_parent_particle : data)
						{
							if (candidate_parent_particle.unique_tag == parent_unique_tag)
							{
								current_particle = &candidate_parent_particle;
								break;
							}
						}

						if (current_particle == nullptr)
							throw std::runtime_error("Bug! Looking for a parent particle that has been deleted!");
					}	
				}

				// Phew. Send data to stdout.
				// It is possible that this particle was detected without any events,
				// in which case deepest_energy is still -1. We should not print anything in that case
				if (deepest_energy > 0)
				{
					std::cout << deepest_z << '\t'
						<< deepest_energy << '\t'
						<< this_particle.particle_data.kin_energy << '\n';
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

};

}} // namespace nbl::drivers

#endif // __CASCADE_SAVING_PARTICLE_MANAGER_H_
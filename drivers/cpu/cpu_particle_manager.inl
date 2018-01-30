#include <algorithm>

namespace nbl { namespace drivers {

template<typename material_manager_t, typename additional_data>
cpu_particle_manager<material_manager_t, additional_data>
	cpu_particle_manager<material_manager_t, additional_data>::create()
{
	cpu_particle_manager<material_manager_t, additional_data> manager{};
	return manager;
}

template<typename material_manager_t, typename additional_data>
void cpu_particle_manager<material_manager_t, additional_data>::destroy(
	cpu_particle_manager<material_manager_t, additional_data> & manager)
{
	manager.data.clear();
	manager.data.shrink_to_fit();
}

template<typename material_manager_t, typename additional_data>
auto cpu_particle_manager<material_manager_t, additional_data>::push(
	particle* particles, primary_tag_t* tags, particle_index_t N)
-> particle_index_t
{
	data.reserve(data.size() + N);
	for (particle_index_t i = 0; i < N; ++i)
	{
		data.push_back({
			NO_EVENT,
			0,
			-123,  // TODO: vacuum
			particles[i],
			tags[i],
			nullptr
		});
	}
	return N;
}

template<typename material_manager_t, typename additional_data>
template<typename detect_function>
void cpu_particle_manager<material_manager_t, additional_data>::flush_detected(
	detect_function func)
{
	for (auto& this_particle : data)
	{
		if (this_particle.status == DETECTED)
		{
			func(this_particle.particle_data, this_particle.primary_tag);
			this_particle.status = TERMINATED;
		}
	}
}

template<typename material_manager_t, typename additional_data>
void cpu_particle_manager<material_manager_t, additional_data>::flush_terminated()
{
	data.erase
	(
		std::remove_if(data.begin(), data.end(),
			[](particle_struct<additional_data> const & x) -> bool
			{ return x.status == TERMINATED; }),
		data.end()
	);
}

template<typename material_manager_t, typename additional_data>
auto cpu_particle_manager<material_manager_t, additional_data>::get_total_count() const -> particle_index_t
{
	return data.size();
}
template<typename material_manager_t, typename additional_data>
auto cpu_particle_manager<material_manager_t, additional_data>::get_running_count() const -> particle_index_t
{
	return static_cast<particle_index_t>(
	std::count_if(data.begin(), data.end(),
		[](particle_struct<additional_data> const & x) -> bool
		{ return x.status != TERMINATED && x.status != DETECTED; }));
}
template<typename material_manager_t, typename additional_data>
auto cpu_particle_manager<material_manager_t, additional_data>::get_detected_count() const -> particle_index_t
{
	return std::count_if(data.begin(), data.end(),
		[](particle_struct<additional_data> const & x) -> bool
		{ return x.status == DETECTED; });
}

template<typename material_manager_t, typename additional_data>
PHYSICS particle & cpu_particle_manager<material_manager_t, additional_data>::operator[](particle_index_t i)
{
	return data[i].particle_data;
}
template<typename material_manager_t, typename additional_data>
PHYSICS particle const & cpu_particle_manager<material_manager_t, additional_data>::operator[](particle_index_t i) const
{
	return data[i].particle_data;
}

template<typename material_manager_t, typename additional_data>
PHYSICS bool cpu_particle_manager<material_manager_t, additional_data>::exists(particle_index_t i) const
{
	return i < data.size();
}

template<typename material_manager_t, typename additional_data>
PHYSICS bool cpu_particle_manager<material_manager_t, additional_data>::active(
	particle_index_t i) const
{
	switch (data[i].status)
	{
		case TERMINATED:
		case DETECTED:
			return false;
		default:
			return true;
	}
}

template<typename material_manager_t, typename additional_data>
PHYSICS auto cpu_particle_manager<material_manager_t, additional_data>::get_material_index(particle_index_t i) const
-> material_index_t
{
	return data[i].current_material;
}
template<typename material_manager_t, typename additional_data>
PHYSICS void cpu_particle_manager<material_manager_t, additional_data>::set_material_index(
	particle_index_t particle_idx, material_index_t new_material_idx)
{
	data[particle_idx].current_material = new_material_idx;
}

template<typename material_manager_t, typename additional_data>
PHYSICS triangle const * cpu_particle_manager<material_manager_t, additional_data>::get_last_triangle(
	particle_index_t i) const
{
	return data[i].last_triangle;
}
template<typename material_manager_t, typename additional_data>
PHYSICS void cpu_particle_manager<material_manager_t, additional_data>::forget_last_triangle(
	particle_index_t i)
{
	data[i].last_triangle = nullptr;
}

template<typename material_manager_t, typename additional_data>
PHYSICS bool cpu_particle_manager<material_manager_t, additional_data>::next_scatter(
	particle_index_t i) const
{
	return data[i].status == SCATTER_EVENT;
}
template<typename material_manager_t, typename additional_data>
PHYSICS uint8_t cpu_particle_manager<material_manager_t, additional_data>::get_next_scatter(
	particle_index_t i) const
{
	return data[i].next_scatter;
}
template<typename material_manager_t, typename additional_data>
PHYSICS bool cpu_particle_manager<material_manager_t, additional_data>::next_intersect(particle_index_t i) const
{
	return data[i].status == INTERSECT_EVENT;
}

template<typename material_manager_t, typename additional_data>
PHYSICS void cpu_particle_manager<material_manager_t, additional_data>::create_secondary(
	particle_index_t primary_idx, particle secondary_particle)
{
	data.push_back({
		NEW_SECONDARY,
		0,
		get_material_index(primary_idx),
		secondary_particle,
		data[primary_idx].primary_tag,
		nullptr
	});
}

template<typename material_manager_t, typename additional_data>
PHYSICS void cpu_particle_manager<material_manager_t, additional_data>::terminate(particle_index_t i)
{
	data[i].status = TERMINATED;
}
template<typename material_manager_t, typename additional_data>
PHYSICS void cpu_particle_manager<material_manager_t, additional_data>::detect(particle_index_t i)
{
	data[i].status = DETECTED;
}

// TODO: the two functions below recalculate the normalization of "dir"
// which has already been done by the driver...
template<typename material_manager_t, typename additional_data>
PHYSICS void cpu_particle_manager<material_manager_t, additional_data>::set_scatter_event(
	particle_index_t i, scatter_event event)
{
	if (event.type != 0)
	{
		data[i].status = SCATTER_EVENT;
		data[i].next_scatter = event.type;
	}
	else
	{
		data[i].status = NO_EVENT;
	}
	data[i].particle_data.pos += normalised(data[i].particle_data.dir) * event.distance;
}
template<typename material_manager_t, typename additional_data>
PHYSICS void cpu_particle_manager<material_manager_t, additional_data>::set_intersect_event(
	particle_index_t i, intersect_event event)
{
	data[i].status = INTERSECT_EVENT;
	data[i].last_triangle = event.isect_triangle;
	data[i].particle_data.pos += normalised(data[i].particle_data.dir) * event.isect_distance;
}

}} // namespace nbl::drivers

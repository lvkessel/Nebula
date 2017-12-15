namespace nbl { namespace drivers {

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::cpu_driver(
	geometry_manager_t geometry,
	intersect_t intersect,
	std::vector<material_t> materials
) :
	_particles(particle_manager_t::create()),
	_materials(material_manager_t::create(materials)),
	_geometry(geometry),
	_intersect(intersect)
{
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::~cpu_driver()
{
	particle_manager_t::destroy(_particles);
	material_manager_t::destroy(_materials);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::push(
	particle* particles,
	primary_tag_t* tags,
	particle_index_t N
) -> particle_index_t
{
	return _particles.push(particles, tags, N);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::do_iteration()
{
	const auto particle_count = _particles.get_total_count();
	for (particle_index_t particle_idx = 0; particle_idx < particle_count; ++particle_idx)
	{
		if (!_particles.active(particle_idx))
			continue;
		
		init(particle_idx);
		intersect(particle_idx);
		scatter(particle_idx);
	}
	_particles.flush_terminated();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::simulate_until_end()
{
	for (particle_index_t particle_idx = 0; particle_idx < _particles.get_total_count(); ++particle_idx)
	{
		while (_particles.active(particle_idx))
		{
			init(particle_idx);
			intersect(particle_idx);
			scatter(particle_idx);
		}
	}
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::get_running_count() const
-> particle_index_t
{
	return _particles.get_running_count();
}
template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
auto cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::get_detected_count() const
-> particle_index_t
{
	return _particles.get_detected_count();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
template<typename detect_function>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::flush_detected(
	detect_function func)
{
	_particles.flush_detected(func);
}
template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::flush_terminated()
{
	_particles.flush_terminated();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::init(particle_index_t particle_idx)
{
	// Get data from memory
	auto this_particle = _particles[particle_idx];

	// If not in domain, terminate
	if (!_geometry.in_domain(this_particle.pos))
	{
		_particles.terminate(particle_idx);
		return;
	}

	// Next scattering event
	scatter_event next_scatter{
		0,
		_geometry.get_max_extent()
	};

	// If not in a vacuum, get next scatter event
	auto this_material_idx = _particles.get_material_index(particle_idx);
	if (_materials.is_physical(this_material_idx))
	{
		const auto this_material = _materials[this_material_idx];

		// Terminate if we can't reach the vacuum
		// It is more natural to put this further up.
		if (!this_material.can_reach_vacuum(this_particle.kin_energy))
		{
			_particles.terminate(particle_idx);
			return;
		}

		// Sample next scattering event
		// TODO: case of no scattering events!
		next_scatter = this_material.sample_path(this_particle, rand_state);
	}

	// Move particle to next event, unless there is a triangle in the way
	this_particle.dir = normalised(this_particle.dir);
	intersect_event next_intersect = _geometry.propagate(
		this_particle.pos, this_particle.dir, next_scatter.distance,
		_particles.get_last_triangle(particle_idx),
		_particles.get_material_index(particle_idx)
	);
	
	if (next_intersect.isect_triangle == nullptr)
	{
		// No triangle intersection: move to scattering position.
		// Scatter there later (after sorting)
		_particles.set_scatter_event(particle_idx, next_scatter);
	}
	else
	{
		// Triangle intersection: move to triangle position.
		_particles.set_intersect_event(particle_idx, next_intersect);
	}
}


template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::intersect(particle_index_t particle_idx)
{
	// ignore all particles except those with an intersect event.
	if (!_particles.next_intersect(particle_idx))
		return;

	_intersect.execute(_materials, _particles, particle_idx, rand_state);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t,
	template<typename> class particle_manager
>
void cpu_driver<scatter_list_t, intersect_t, geometry_manager_t, particle_manager>::scatter(particle_index_t particle_idx)
{
	// ignore all particles except those with an inelastic event.
	if (!_particles.next_scatter(particle_idx))
		return;

	// forget last intersected triangle. This event might cause us to scatter back into that triangle
	// and we don't want to ignore that triangle if so.
	_particles.forget_last_triangle(particle_idx);

	_materials[_particles.get_material_index(particle_idx)].execute(
		_particles.get_next_scatter(particle_idx), _particles, particle_idx, rand_state);
}

}} // namespace nbl::drivers

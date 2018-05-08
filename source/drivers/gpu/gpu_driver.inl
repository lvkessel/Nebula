namespace nbl { namespace drivers {

namespace kernels
{
	__global__ void init_random_states(util::random_generator<true>* curand_states, unsigned long long seed, size_t capacity);

	template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t>
	__global__ void init(particle_manager_t particles,
		material_manager_t materials,
		geometry_manager_t geometry,
		util::random_generator<true>* curand_states);

	template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t, typename intersect_t>
	__global__ void intersect(particle_manager_t particles,
		material_manager_t materials,
		geometry_manager_t geometry,
		util::random_generator<true>* curand_states,
		intersect_t isect);

	template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t>
	__global__ void inelastic(particle_manager_t particles,
		material_manager_t materials,
		geometry_manager_t geometry,
		util::random_generator<true>* curand_states);

	template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t>
	__global__ void elastic(particle_manager_t particles,
		material_manager_t materials,
		geometry_manager_t geometry,
		util::random_generator<true>* curand_states);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::gpu_driver(
	particle_index_t particle_capacity,
	geometry_manager_t geometry,
	intersect_t intersect,
	std::vector<material_t> materials,
	seed_t seed
) :
	_particles(particle_manager_t::create(particle_capacity)),
	_materials(material_manager_t::create(materials)),
	_geometry(geometry),
	_intersect(intersect),
	_num_blocks(1 + particle_capacity/_threads_per_block)
{
	/*
	 * Init random states
	 */
	cuda::cuda_new<util::random_generator<true>>(&curand_states, particle_capacity);
	kernels::init_random_states<<<_num_blocks, _threads_per_block>>>(
		curand_states, seed, particle_capacity
	);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::~gpu_driver()
{
	particle_manager_t::destroy(_particles);
	material_manager_t::destroy(_materials);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU auto gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::push(
	particle* particles,
	uint32_t* tags,
	particle_index_t N
) -> particle_index_t
{
	return _particles.push(particles, tags, N);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU void gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::do_iteration()
{
	init();
	_particles.sort();
	events();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU auto gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::get_running_count() const
-> particle_index_t
{
	return _particles.get_running_count();
}
template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU auto gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::get_detected_count() const
-> particle_index_t
{
	return _particles.get_detected_count();
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
template<typename detect_function>
CPU void gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::flush_detected(
	detect_function func)
{
	_particles.flush_detected(func);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU void gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::init()
{
	kernels::init<<<_num_blocks, _threads_per_block>>>(
		_particles, _materials, _geometry, curand_states
	);
}

template<typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t
>
CPU void gpu_driver<scatter_list_t, intersect_t, geometry_manager_t>::events()
{
	kernels::intersect<<<_num_blocks, _threads_per_block>>>(
		_particles, _materials, _geometry, curand_states, _intersect
	);
	kernels::elastic<<<_num_blocks, _threads_per_block>>>(
		_particles, _materials, _geometry, curand_states
	);
	kernels::inelastic<<<_num_blocks, _threads_per_block>>>(
		_particles, _materials, _geometry, curand_states
	);
}

__global__ void kernels::init_random_states(
	util::random_generator<true>* curand_states,
	unsigned long long seed,
	size_t capacity)
{
	const auto i = threadIdx.x+blockIdx.x*blockDim.x;
	if(i >= capacity)
		return;

	curand_states[i] = util::random_generator<true>(seed, i);
}

template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t>
__global__ void kernels::init(
	particle_manager_t particles,
	material_manager_t materials,
	geometry_manager_t geometry,
	util::random_generator<true>* curand_states)
{
	const auto particle_idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(!particles.exists(particle_idx))
		return;

	// Ignore all particles that are not active
	if (!particles.active(particle_idx))
		return;

	// Get data from memory
	auto this_particle = particles[particle_idx];

	// If not in domain, terminate
	if (!geometry.in_domain(this_particle.pos))
	{
		particles.terminate(particle_idx);
		return;
	}

	// Next scattering event
	scatter_event next_scatter{
		0,
		geometry.get_max_extent()
	};

	// If not in a vacuum, get next scatter event
	auto this_material_idx = particles.get_material_index(particle_idx);
	if (materials.is_physical(this_material_idx))
	{
		const auto this_material = materials[this_material_idx];

		// Terminate if we can't reach the vacuum
		// It is more natural to put this further up.
		if (!this_material.can_reach_vacuum(this_particle.kin_energy))
		{
			particles.terminate(particle_idx);
			return;
		}

		// Sample next scattering event
		// TODO: case of no scattering events!
		auto rng = curand_states[particle_idx];
		next_scatter = this_material.sample_path(this_particle, rng);
		curand_states[particle_idx] = rng;
	}

	// Move particle to next event, unless there is a triangle in the way
	this_particle.dir = normalised(this_particle.dir);
	intersect_event next_intersect = geometry.propagate(
		this_particle.pos, this_particle.dir, next_scatter.distance,
		particles.get_last_triangle(particle_idx),
		particles.get_material_index(particle_idx)
	);

	if (next_intersect.isect_triangle == nullptr)
	{
		// No triangle intersection: move to scattering position.
		// Scatter there later (after sorting)
		this_particle.pos += this_particle.dir * next_scatter.distance;
		particles.set_scatter_event(particle_idx, next_scatter.type);
	}
	else
	{
		// Triangle intersection: move to triangle position.
		this_particle.pos += this_particle.dir*next_intersect.isect_distance;
		particles.set_intersect_event(particle_idx, next_intersect.isect_triangle);
	}

	// Store new particle data in global memory
	particles[particle_idx] = this_particle;
}

template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t, typename intersect_t>
__global__ void kernels::intersect(
	particle_manager_t particles,
	material_manager_t materials,
	geometry_manager_t geometry,
	util::random_generator<true>* curand_states,
	intersect_t isect)
{
	// Get particle index
	const auto i = threadIdx.x+blockIdx.x*blockDim.x;
	if(!particles.exists(i))
		return;
	const auto particle_idx = particles.get_particle_index(i);

	// ignore all particles except those with an intersect event.
	if (!particles.next_intersect(particle_idx))
		return;

	auto rng = curand_states[particle_idx];
	isect.execute(materials, particles, particle_idx, rng);
	curand_states[particle_idx] = rng;
}

template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t>
__global__ void kernels::inelastic(
	particle_manager_t particles,
	material_manager_t materials,
	geometry_manager_t geometry,
	util::random_generator<true>* curand_states)
{
	// Get particle index
	const auto i = threadIdx.x+blockIdx.x*blockDim.x;
	if(!particles.exists(i))
		return;
	const auto particle_idx = particles.get_particle_index(i);

	// ignore all particles except those with an inelastic event.
	if (!particles.next_inelastic(particle_idx))
		return;

	// If we can't make a secondary right now, set to pending and wait an iteration.
	if (!particles.secondary_slot_free())
	{
		particles.pending(particle_idx);
		return;
	}

	// We are definitely scattering inelastically now (status needs to be set if previously pending)
	particles.inelastic(particle_idx);

	// forget last intersected triangle. This event might cause us to scatter back into that triangle
	// and we don't want to ignore that triangle if so.
	particles.forget_last_triangle(particle_idx);

	// TODO: event types hardcoded here
	auto rng = curand_states[particle_idx];
	materials[particles.get_material_index(particle_idx)].execute(1, particles, particle_idx, rng);
	curand_states[particle_idx] = rng;
}

template<typename particle_manager_t, typename material_manager_t, typename geometry_manager_t>
__global__ void kernels::elastic(particle_manager_t particles,
	material_manager_t materials,
	geometry_manager_t geometry,
	util::random_generator<true>* curand_states)
{
	// Get particle index
	const auto i = threadIdx.x+blockIdx.x*blockDim.x;
	if(!particles.exists(i))
		return;
	const auto particle_idx = particles.get_particle_index(i);

	// ignore all particles except those with an elastic event.
	if (!particles.next_elastic(particle_idx))
		return;

	// forget last intersected triangle. This event might cause us to scatter back into that triangle
	// and we don't want to ignore that triangle if so.
	particles.forget_last_triangle(particle_idx);

	// TODO: event types hardcoded here
	auto rng = curand_states[particle_idx];
	materials[particles.get_material_index(particle_idx)].execute(2, particles, particle_idx, rng);
	curand_states[particle_idx] = rng;
}

}} // namespace nbl::drivers

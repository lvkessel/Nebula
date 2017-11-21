#include "config/config.h"
#include "core/material.h"
#include "physics/inelastic_thomas.h"
#include "physics/elastic_thomas.h"
#include "physics/intersect_thomas.h"

#include "drivers/gpu/gpu_driver.h"
#include "drivers/cpu/cpu_driver.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include "legacy_thomas/load_tri_file.h"
#include "legacy_thomas/load_pri_file.h"
#include "legacy_thomas/load_mat_file.h"

#define USE_GPU true

// Main typedefs
using geometry_t = nbl::geometry::octree<USE_GPU>;

using thomas_scatter_list = scatter_list<
	nbl::scatter::inelastic_thomas<USE_GPU>,
	nbl::scatter::elastic_thomas<USE_GPU>
>;
using thomas_material = material<thomas_scatter_list>;
using intersect_t = intersect_thomas<>;

#if USE_GPU
using driver = nbl::drivers::gpu_driver<
	thomas_scatter_list,
	intersect_t,
	geometry_t
>;
#else
using driver = nbl::drivers::cpu_driver<
	thomas_scatter_list,
	intersect_t,
	geometry_t
>;
#endif

// TODO: material not really destroyed.
thomas_material load_material(std::string const & filename)
{
	material_legacy_thomas mat_legacy = load_mat_file(filename);
	return thomas_material(mat_legacy);
}

int main(int argc, char** argv)
{
	// Settings
	const size_t capacity = 200000;
	const size_t prescan_size = 1000;
	const real batch_factor = .9_r;

	const std::string usage(std::string("Usage: ") + argv[0] +
		"<geometry.tri> <primaries.pri> [material0.mat] .. [materialN.mat]");

	if (argc <= 3)
	{
		std::clog << usage << std::endl;
		return 1;
	}

	// Load geometry
	std::vector<triangle> triangles = load_tri_file(argv[1]);
	if (triangles.empty())
	{
		std::clog << "Error: could not load triangles!\n" << usage << std::endl;
		return 1;
	}
	geometry_t geometry = geometry_t::create(triangles);


	// Load primaries
	std::vector<particle> primaries; std::vector<int2> pixels;
	std::tie(primaries, pixels) = separate_pairs(
		sort_pri_file(load_pri_file(argv[2]), prescan_size));
	if (primaries.empty())
	{
		std::clog << "Error: could not load primary electrons!\n" << usage << std::endl;
		return 1;
	}

	// The GPU driver only accepts uint32 tags. So we make a map: GPU tag is
	// the index of the primary particle in the "primaries" / "pixels" array.
	std::vector<uint32_t> gpu_tags(primaries.size());
	std::iota(gpu_tags.begin(), gpu_tags.end(), 0); // Fill with 0, 1, ... tags.size()-1

	// These variables are used for the exposure
	particle* next_primary = primaries.data();
	uint32_t* next_tag = gpu_tags.data();
	size_t primaries_to_go = primaries.size();

	// Load materials
	std::vector<thomas_material> materials;
	for (int parameter_idx = 3; parameter_idx < argc; ++parameter_idx)
		materials.push_back(load_material(argv[parameter_idx]));

	intersect_t inter;

	// Print debug data
	std::clog << "Loaded " << triangles.size() << " triangles." << std::endl
		<< "  min = {" << geometry.AABB_min().x << ", " << geometry.AABB_min().y << ", " << geometry.AABB_min().z << "}" << std::endl
		<< "  max = {" << geometry.AABB_max().x << ", " << geometry.AABB_max().y << ", " << geometry.AABB_max().z << "}" << std::endl;
	std::clog << "Loaded " << primaries.size() << " primaries." << std::endl;
	std::clog << "Loaded " << materials.size() << " materials." << std::endl;

	// Prepare driver
#if USE_GPU
	driver d(capacity, geometry, inter, materials);
#else
	driver d(geometry, inter, materials);
#endif

	// First, do the prescan
	std::vector<std::pair<uint32_t, uint32_t>> prescan_stats; // Holds (running_count, detected_count)
	// Push first batch
	{
		auto particles_pushed = d.push(next_primary, next_tag, std::min(prescan_size, primaries_to_go));
		next_primary += particles_pushed;
		next_tag += particles_pushed;
		primaries_to_go -= particles_pushed;
		prescan_stats.push_back({ particles_pushed, 0 });
	}
	// Execute prescan
	while (prescan_stats.back().first > 0)
	{
		d.do_iteration();

		// TODO: this can be optimised to just one function with one loop.
		prescan_stats.push_back({ d.get_running_count(), d.get_detected_count() });
		std::clog << " \rExecuting pre-scan"
			<< " | running: " << prescan_stats.back().first
			<< " | detected: " << prescan_stats.back().second;
	}
	// Find frame_size and batch_size based on the prescan stats.
	// frame_size is the iteration number where running_count was maximal.
	const size_t frame_size = 1 + std::distance(prescan_stats.begin(),
		std::max_element(prescan_stats.begin(), prescan_stats.end(),
		[](std::pair<uint32_t, uint32_t> p1, std::pair<uint32_t, uint32_t> p2) -> bool
		{ return p1.first < p2.first; }));
	// Batch size
	size_t batch_size;
	{
		real accumulator = 0;
		for (size_t i = 2*frame_size; i < prescan_stats.size(); i += frame_size)
			accumulator += prescan_stats[i].first / real(prescan_size);
		accumulator += 2*prescan_stats[frame_size].first / real(prescan_size);
		accumulator += 2*prescan_stats[frame_size].second / real(prescan_size);
		batch_size = size_t(batch_factor*capacity / accumulator);
	}
	std::clog << "\nframe_size = " << frame_size << " | batch_size = " << batch_size << std::endl;

	std::ofstream of("tmp.bin", std::ofstream::binary);

	for (;;)
	{
		// Push new batch
		{
			auto particles_pushed = d.push(next_primary, next_tag, std::min(batch_size, primaries_to_go));
			next_primary += particles_pushed;
			next_tag += particles_pushed;
			primaries_to_go -= particles_pushed;
		}

		// Execute frame
		for (uint32_t i = 0; i < frame_size; ++i)
			d.do_iteration();
#if USE_GPU
		cudaDeviceSynchronize();
#endif

		// Flush output data
		d.flush_detected([&of, &pixels](particle p, uint32_t t)
		{
			float buffer[7];
			buffer[0] = p.pos.x; buffer[1] = p.pos.y; buffer[2] = p.pos.z;
			buffer[3] = p.dir.x; buffer[4] = p.dir.y; buffer[5] = p.dir.z;
			buffer[6] = p.kin_energy;
			of.write(reinterpret_cast<const char*>(buffer), sizeof(buffer));
			of.write(reinterpret_cast<const char*>(&pixels[t]), sizeof(int2));
/*
			// Write out
			of.write(reinterpret_cast<const char*>(&(p.pos.x)), sizeof(p.pos.x));
			of.write(reinterpret_cast<const char*>(&(p.pos.y)), sizeof(p.pos.y));
			of.write(reinterpret_cast<const char*>(&(p.pos.z)), sizeof(p.pos.z));
			of.write(reinterpret_cast<const char*>(&(p.dir.x)), sizeof(p.dir.x));
			of.write(reinterpret_cast<const char*>(&(p.dir.y)), sizeof(p.dir.y));
			of.write(reinterpret_cast<const char*>(&(p.dir.z)), sizeof(p.dir.z));
			of.write(reinterpret_cast<const char*>(&(p.kin_energy)), sizeof(p.kin_energy));
*/
		});

		// Show progress, end when finished.
		auto running_count = d.get_running_count();

		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - ((double)primaries_to_go / primaries.size()))
			<< "%, # Running: " << running_count;
		if (running_count == 0 && primaries_to_go == 0)
			break;
	}

#if USE_GPU
	cudaError_t err = cudaDeviceSynchronize();
	std::clog << std::endl << "CUDA error code " << err << std::endl;
#endif

	std::clog << "Done." << std::endl;
	//std::cin.get();
	return 0;
}
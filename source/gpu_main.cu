#include "config/config.h"
#include "physics_config.h"
#include "core/material.h"
#include "common/cli_params.h"

#include "drivers/gpu/gpu_driver.h"
#include "drivers/cpu/cpu_driver.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <numeric>
#include "legacy_thomas/load_tri_file.h"
#include "legacy_thomas/load_pri_file.h"
#include "legacy_thomas/load_mat_file.h"

// Main typedefs
using geometry_t = nbl::geometry::octree<true>;
using material_t = material<scatter_physics<true>>;

using driver = nbl::drivers::gpu_driver<
	scatter_physics<true>,
	intersect_t,
	geometry_t
>;

// TODO: material not really destroyed.
material_t load_material(std::string const & filename)
{
	if (filename.back() == 't')
	{
		// Old .mat file format
		return material_t(load_mat_file(filename));
	}
	else
	{
		// New HDF5 file format
		return material_t(nbl::hdf5_file(filename));
	}
}

struct timelog
{
	void add(std::string const & name,
		std::chrono::steady_clock::time_point t1,
		std::chrono::steady_clock::time_point t2)
	{
		data.push_back({name, t2-t1});
	}

	void print(std::ostream& out)
	{
		const size_t len_col1 = std::max_element(data.begin(), data.end(),
			[](data_t const & d1, data_t const & d2) -> bool
			{ return d1.first.size() < d2.first.size(); })->first.size();

		out << "Timing (seconds):\n";
		for (auto&& d : data)
		{
			const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(d.second);
			out << std::setw(len_col1) << d.first << ": "
				<< std::fixed << std::setprecision(3) << (ms.count()/1000.)
				<< '\n';
		}
	}

	using data_t = std::pair<std::string, std::chrono::steady_clock::duration>;
	std::vector<data_t> data;
};

int main(int argc, char** argv)
{
	// Settings
	size_t capacity = 1000000;
	size_t prescan_size = 1000;
	real batch_factor = .9_r;
	typename driver::seed_t seed = 0x14f8214e78c7e39b;
	bool sort_primaries = false;

	cli_params p(argc, argv);
	p.get_optional_flag("capacity", capacity);
	p.get_optional_flag("prescan-size", prescan_size);
	p.get_optional_flag("batch-factor", batch_factor);
	p.get_optional_flag("seed", seed);
	p.get_optional_flag("sort-primaries", sort_primaries);

	const std::string usage("Usage: " + p.get_program_name() +
		" [options] <geometry.tri> <primaries.pri> [material0.mat] .. [materialN.mat]\n"
		"Options:\n"
		"\t--capacity       [1000000]\n"
		"\t--prescan-size   [1000]\n"
		"\t--batch-factor   [0.9]\n"
		"\t--seed           [0x14f8214e78c7e39b]\n"
		"\t--sort-primaries [0]\n");

	// Setup time logging
	timelog timer;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;

	// Interpret command-line options
	std::vector<std::string> pos_flags = p.get_positional();
	if (pos_flags.size() < 3 || capacity <= 0 || prescan_size <= 0 || batch_factor <= 0)
	{
		std::clog << usage << std::endl;
		return 1;
	}

	// Load geometry
	std::clog << "Loading geometry..." << std::endl;
	t1 = std::chrono::steady_clock::now();
	std::vector<triangle> triangles = load_tri_file(pos_flags[0]);
	t2 = std::chrono::steady_clock::now();
	timer.add("Loading triangles", t1, t2);

	if (triangles.empty())
	{
		std::clog << "Error: could not load triangles!\n" << usage << std::endl;
		return 1;
	}

	t1 = std::chrono::steady_clock::now();
	geometry_t geometry = geometry_t::create(triangles);
	t2 = std::chrono::steady_clock::now();
	timer.add("Building acceleration structure", t1, t2);


	// Load primaries
	std::clog << "Loading primary electrons..." << std::endl;
	t1 = std::chrono::steady_clock::now();
	std::vector<particle> primaries; std::vector<int2> pixels;
	std::tie(primaries, pixels) = load_pri_file(pos_flags[1], geometry.AABB_min(), geometry.AABB_max());
	if (sort_primaries)
		sort_pri_file(primaries, pixels);
	prescan_shuffle(primaries, pixels, prescan_size);
	t2 = std::chrono::steady_clock::now();
	timer.add("Loading primary electrons", t1, t2);

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
	std::clog << "Loading materials..." << std::endl;
	t1 = std::chrono::steady_clock::now();
	std::vector<material_t> materials;
	for (size_t parameter_idx = 2; parameter_idx < pos_flags.size(); ++parameter_idx)
		materials.push_back(load_material(pos_flags[parameter_idx]));
	t2 = std::chrono::steady_clock::now();
	timer.add("Loading materials", t1, t2);


	intersect_t inter;

	// Print debug data
	std::clog << "\n"
		<< "Loaded " << triangles.size() << " triangles.\n"
		<< "  min = {" << geometry.AABB_min().x << ", " << geometry.AABB_min().y << ", " << geometry.AABB_min().z << "}\n"
		<< "  max = {" << geometry.AABB_max().x << ", " << geometry.AABB_max().y << ", " << geometry.AABB_max().z << "}\n"
		<< "Loaded " << primaries.size() << " primaries.\n"
		<< "Loaded " << materials.size() << " materials.\n\n" << std::flush;

	// Prepare driver
	driver d(capacity, geometry, inter, materials, seed);

	// First, do the prescan
	t1 = std::chrono::steady_clock::now();
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
	t2 = std::chrono::steady_clock::now();
	timer.add("Prescan", t1, t2);


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
	d.allocate_input_buffers(batch_size);

	//std::ofstream of("tmp.bin", std::ofstream::binary);
	std::ostream& of = std::cout;

	// Simulation
	t1 = std::chrono::steady_clock::now();
	for (;;)
	{
		// Copy the simulation to buffer,
		// push new data into the simulation
		d.buffer_detected();
		d.push_to_simulation();
		cudaDeviceSynchronize();

		// Execute frame
		for (uint32_t i = 0; i < frame_size; ++i)
			d.do_iteration();

		// Push new batch asynchronously
		{
			auto particles_pushed = d.push_to_buffer(next_primary, next_tag, std::min(batch_size, primaries_to_go));
			next_primary += particles_pushed;
			next_tag += particles_pushed;
			primaries_to_go -= particles_pushed;
		}

		// Output detected electrons from buffer
		auto running_count = d.flush_buffered([&of, &pixels](particle p, uint32_t t)
		{
			float buffer[7];
			buffer[0] = p.pos.x; buffer[1] = p.pos.y; buffer[2] = p.pos.z;
			buffer[3] = p.dir.x; buffer[4] = p.dir.y; buffer[5] = p.dir.z;
			buffer[6] = p.kin_energy;
			of.write(reinterpret_cast<const char*>(buffer), sizeof(buffer));
			of.write(reinterpret_cast<const char*>(&pixels[t]), sizeof(int2));
		});

		// Show progress
		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - ((double)primaries_to_go / primaries.size()))
			<< "%, # Running: " << running_count;

		cudaDeviceSynchronize();

		if (running_count == 0 && primaries_to_go == 0)
			break;
	}
	t2 = std::chrono::steady_clock::now();
	timer.add("Simulation", t1, t2);

	cudaError_t err = cudaDeviceSynchronize();
	if (err == 0)
		std::clog << "\nSimulation successful!\n\n";
	else
		std::clog << "\nSimulation ended with CUDA error code " << err << "\n\n";

	timer.print(std::clog);
	return 0;
}

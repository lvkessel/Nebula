#include "config/config.h"
#include "core/material.h"
#include "physics/inelastic_thomas.h"
#include "physics/inelastic_penn.h"
#include "physics/elastic_thomas.h"
#include "physics/intersect_thomas.h"

#include "drivers/cpu/cpu_driver.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "legacy_thomas/load_tri_file.h"
#include "legacy_thomas/load_pri_file.h"
#include "legacy_thomas/load_mat_file.h"

#include <thread>
#include <mutex>

#define TIMING 1

// Main typedefs
using geometry_t = nbl::geometry::octree<false>;

using scatter_physics = scatter_list<
	//nbl::scatter::inelastic_thomas<false>,
	nbl::scatter::inelastic_penn<false>,
	nbl::scatter::elastic_thomas<false>
>;
using material_t = material<scatter_physics>;
using intersect_t = intersect_thomas<>;

using driver = nbl::drivers::cpu_driver<
	scatter_physics,
	intersect_t,
	geometry_t>;

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

class work_pool
{
public:
	work_pool(particle* primaries, uint32_t* tags, size_t N) :
		next_primary(primaries), next_tag(tags), primaries_to_go(N)
	{}

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

	size_t get_primaries_to_go()
	{
		std::lock_guard<std::mutex> lock(mutex);
		return primaries_to_go;
	}

private:
	std::mutex mutex;

	particle* next_primary;
	uint32_t* next_tag;
	size_t primaries_to_go;
};

class output_stream
{
public:
	output_stream(std::string const & filename, std::vector<int2> pixel_map)
		: file(filename, std::ofstream::binary), pixel_map(pixel_map)
	{}

	void flush(driver& d)
	{
		std::lock_guard<std::mutex> lock(mutex);
		d.flush_detected([this](particle p, uint32_t t)
		{
			// Simulation might have been done in single or double precision.
			// e-scatter's output files are always single precision, so copy to buffer like this.
			float buffer[7];
			buffer[0] = p.pos.x; buffer[1] = p.pos.y; buffer[2] = p.pos.z;
			buffer[3] = p.dir.x; buffer[4] = p.dir.y; buffer[5] = p.dir.z;
			buffer[6] = p.kin_energy;
			file.write(reinterpret_cast<const char*>(buffer), sizeof(buffer));
			file.write(reinterpret_cast<const char*>(&(pixel_map[t])), sizeof(int2));
		});
	}

private:
	std::ofstream file;
	std::mutex mutex;
	std::vector<int2> pixel_map;
};

int main(int argc, char** argv)
{
	// Settings
	// TODO: all of this is completely pointless for CPU.
	const size_t capacity = 200000;
	const size_t prescan_size = 1000;
	const real batch_factor = .9_r;

	const std::string usage(std::string("Usage: ") + argv[0] + "<geometry.tri> <primaries.pri> [material0.mat] .. [materialN.mat]");

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

	// This manages the work to be done (thread-safe).
	work_pool pool(primaries.data(), gpu_tags.data(), primaries.size());
	// And this represents the output stream (thread-safe)
	output_stream out_file("tmp.bin", pixels);


	// Load materials
	std::vector<material_t> materials;
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
	driver d0(geometry, inter, materials);

	// First, do the prescan
	std::vector<std::pair<uint32_t, uint32_t>> prescan_stats; // Holds (running_count, detected_count)
	// Push first batch
	{
		auto work_data = pool.get_work(prescan_size);
		auto particles_pushed = d0.push(
			std::get<0>(work_data),  // particle*
			std::get<1>(work_data),  // tag*
			std::get<2>(work_data)); // number
		prescan_stats.push_back({ particles_pushed, 0 });
	}
	// Execute prescan
	while (prescan_stats.back().first > 0)
	{
		d0.do_iteration();

		// TODO: this can be optimised to just one function with one loop.
		prescan_stats.push_back({ d0.get_running_count(), d0.get_detected_count() });
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

	// Simulation loop
	auto sim_loop = [&pool, &out_file, &geometry, &inter, &materials, frame_size, batch_size]() -> void
	{
		driver d(geometry, inter, materials);
		for (;;)
		{
			// Push new particles
			auto work_data = pool.get_work(batch_size);
			auto particles_pushed = d.push(
				std::get<0>(work_data),  // particle*
				std::get<1>(work_data),  // tag*
				std::get<2>(work_data)); // number

			// Simulate a little
			for (uint32_t i = 0; i < frame_size; ++i)
				d.do_iteration();

			// Flush output data
			out_file.flush(d);

			if (std::get<2>(work_data) == 0 && d.get_running_count() == 0)
				return;
		}
	};

	const auto n_threads = std::thread::hardware_concurrency();
	std::clog << "Creating " << n_threads << " CPU drivers" << std::endl;
	std::vector<std::thread> threads;
#if TIMING
	const auto t1 = std::chrono::high_resolution_clock::now();
#endif
	for (unsigned int i = 0; i < n_threads; ++i)
		threads.push_back(std::thread(sim_loop));

#if !TIMING
	for (;;)
	{
		std::this_thread::sleep_for(std::chrono::seconds(1));
		auto primaries_to_go = pool.get_primaries_to_go();
		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - ((double)primaries_to_go / primaries.size()));
		if (primaries_to_go == 0)
			break;
	}
#endif

	for (auto& t : threads)
		t.join();

#if TIMING
	const auto t2 = std::chrono::high_resolution_clock::now();
	std::clog << (std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.) << std::endl;
#endif

	std::clog << "Done." << std::endl;
	//std::cin.get();
	return 0;
}

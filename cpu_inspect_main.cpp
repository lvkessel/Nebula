#include "config/config.h"
#include "physics_config.h"
#include "core/material.h"

#include "drivers/cpu/cpu_driver.h"
#include "drivers/tagging/cascade_saving_particle_manager.h"

#include "geometry/trilist.h"
#include "geometry/octree.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include "legacy_thomas/load_tri_file.h"
#include "legacy_thomas/load_pri_file.h"
#include "legacy_thomas/load_mat_file.h"

// Main typedefs
using geometry_t = nbl::geometry::octree<false>;
using material_t = material<scatter_physics<false>>;

using driver = nbl::drivers::cpu_driver<
	scatter_physics<false>,
	intersect_t,
	geometry_t,
	nbl::drivers::cascade_saving_particle_manager
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

int main(int argc, char** argv)
{
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
	std::tie(primaries, pixels) = separate_pairs(load_pri_file(argv[2]));
	if (primaries.empty())
	{
		std::clog << "Error: could not load primary electrons!\n" << usage << std::endl;
		return 1;
	}

	// The GPU driver only accepts uint32 tags. So we make a map: GPU tag is
	// the index of the primary particle in the "primaries" / "pixels" array.
	std::vector<uint32_t> gpu_tags(primaries.size());
	std::iota(gpu_tags.begin(), gpu_tags.end(), 0); // Fill with 0, 1, ... tags.size()-1


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
	driver d(geometry, inter, materials);

	std::ofstream of("tmp.bin", std::ofstream::binary);

	// Simulation loop
	for (size_t particle_idx = 0; particle_idx < primaries.size(); ++particle_idx)
	{
		// Push new particle
		auto particles_pushed = d.push(
			primaries.data() + particle_idx, // particle*
			gpu_tags.data() + particle_idx,  // tag*
			1);                              // number

		// Simulate a little
		d.simulate_until_end();

		// Flush output data
		d.flush_detected([&of, &pixels](particle p, uint32_t t)
		{
			float buffer[7];
			buffer[0] = p.pos.x; buffer[1] = p.pos.y; buffer[2] = p.pos.z;
			buffer[3] = p.dir.x; buffer[4] = p.dir.y; buffer[5] = p.dir.z;
			buffer[6] = p.kin_energy;
			of.write(reinterpret_cast<const char*>(buffer), sizeof(buffer));
			of.write(reinterpret_cast<const char*>(&pixels[t]), sizeof(int2));
		});

		d.flush_terminated();

		std::clog << " \rProgress "
			<< std::fixed << std::setprecision(2) << 100 * (1 - (double(primaries.size() - particle_idx) / primaries.size()))
			<< "%";
	}


	std::clog << "Done." << std::endl;
	//std::cin.get();
	return 0;
}

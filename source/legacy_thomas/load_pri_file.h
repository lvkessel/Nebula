#ifndef __LOAD_PRI_FILE_H_
#define __LOAD_PRI_FILE_H_

#include <vector>
#include <tuple>
#include <fstream>
#include <string>
#include "../core/particle.h"

#if !CUDA_HEADERS_AVAILABLE
struct int2 { int x, y; };
#endif

std::vector<std::pair<particle, int2>> load_pri_file(std::string const & filename, vec3 min_pos, vec3 max_pos)
{
	std::vector<std::pair<particle, int2>> particle_vec;

	std::ifstream ifs(filename, std::ifstream::binary);
	if (!ifs.is_open())
		return particle_vec;

	while(!ifs.eof())
	{
		float buffer[7];
		int2 pixel;
		ifs.read(reinterpret_cast<char*>(buffer), sizeof(buffer));
		if (ifs.gcount() == 0)
			break; // end of file
		if (ifs.gcount() != sizeof(buffer))
			throw std::runtime_error("Unexpected end of primaries file");
		ifs.read(reinterpret_cast<char*>(&pixel), sizeof(pixel));
		if (ifs.gcount() != sizeof(pixel))
			throw std::runtime_error("Unexpected end of primaries file");

		particle primary
		{
			{ buffer[0], buffer[1], buffer[2] }, // pos
			{ buffer[3], buffer[4], buffer[5] }, // dir
			buffer[6]                            // energy
		};

		// NAN < xxx == FALSE, so this form works correctly if NaN is read from the file
		if (!(primary.pos.x > min_pos.x && primary.pos.y > min_pos.y && primary.pos.z > min_pos.z &&
			primary.pos.x < max_pos.x && primary.pos.y < max_pos.y && primary.pos.z < max_pos.z))
			std::clog << "WARNING: primary electron starting outside geometry at ("
				<< primary.pos.x << ", " << primary.pos.y << ", " << primary.pos.z << ")." << std::endl;

		if (!(primary.kin_energy > EPSILON))
			std::clog << "WARNING: primary electron starting with low kinetic energy ("
				<< primary.kin_energy << " eV)." << std::endl;

		if (primary.kin_energy > K_max)
			std::clog << "WARNING: primary electron starting with kinetic energy higher than the hard-coded limit ("
				<< (primary.kin_energy/1000) << " keV, max " << (K_max/1000) << " keV). Complain to the developers." << std::endl;

		if (!(std::abs(primary.dir.x) > EPSILON || std::abs(primary.dir.y) > EPSILON || std::abs(primary.dir.z) > EPSILON))
			std::clog << "WARNING: primary electron starts with unphysical direction vector ("
				<< primary.dir.x << ", " << primary.dir.y << ", " << primary.dir.z << ")." << std::endl;

		particle_vec.push_back({ primary, pixel });
	}
	ifs.close(); 
	return particle_vec;
}

uint64_t morton(uint16_t x, uint16_t y, uint16_t z)
{
	uint64_t morton = 0;
	for (unsigned int i = 0; i < 16; i++)
	{
		morton |= (x & (static_cast<uint64_t>(1) << i)) << (2 * i + 0);
		morton |= (y & (static_cast<uint64_t>(1) << i)) << (2 * i + 1);
		morton |= (z & (static_cast<uint64_t>(1) << i)) << (2 * i + 2);
	}
	return morton;
}

/*
 * The first prescan_size particles are random, uniformly shuffled from the exposure
 * for a well-sampled prescan.
 * The other ones are sorted by Morton index, to make sure that nearby particles in
 * the exposure are also spatially close. This gives a significant speedup in the
 * GPU collision detection routine.
 */
std::vector<std::pair<particle, int2>> sort_pri_file(std::vector<std::pair<particle, int2>> data, size_t prescan_size)
{
	if (prescan_size > 0)
		std::shuffle(data.begin(), data.end(), std::default_random_engine());

	if (prescan_size < data.size())
	{
		// Find min, max coordinates of all primaries
		vec3 vmin = data[0].first.pos;
		vec3 vmax = vmin;
		for (const auto& p : data)
		{
			vmin = {
				std::min(vmin.x, p.first.pos.x),
				std::min(vmin.y, p.first.pos.y),
				std::min(vmin.z, p.first.pos.z),
			};
			vmax = {
				std::max(vmax.x, p.first.pos.x),
				std::max(vmax.y, p.first.pos.y),
				std::max(vmax.z, p.first.pos.z),
			};
		}

		// Sort by Morton ordering
		const vec3 vsize = vmax - vmin;
		std::sort(data.begin() + prescan_size, data.end(),
			[vmin, vsize](std::pair<particle, int2> p1, std::pair<particle, int2> p2) -> bool
		{
			const vec3 p1n = (p1.first.pos - vmin);
			const vec3 p2n = (p2.first.pos - vmin);

			const auto m1 = morton(
				static_cast<uint16_t>(p1n.x / vsize.x * std::numeric_limits<uint16_t>::max()),
				static_cast<uint16_t>(p1n.y / vsize.y * std::numeric_limits<uint16_t>::max()),
				static_cast<uint16_t>(p1n.z / vsize.z * std::numeric_limits<uint16_t>::max()));
			const auto m2 = morton(
				static_cast<uint16_t>(p2n.x / vsize.x * std::numeric_limits<uint16_t>::max()),
				static_cast<uint16_t>(p2n.y / vsize.y * std::numeric_limits<uint16_t>::max()),
				static_cast<uint16_t>(p2n.z / vsize.z * std::numeric_limits<uint16_t>::max()));

			return m1 < m2;
		});
	}

	return data;
}

/*
 * Return a pair of vectors (not a vector of pairs), allowing the output
 * of these functions to be fed directly into the simulator.
 */
std::pair<std::vector<particle>, std::vector<int2>> separate_pairs(std::vector<std::pair<particle, int2>> const & data)
{
	std::vector<particle> particle_vec; particle_vec.reserve(data.size());
	std::vector<int2> pixel_vec; pixel_vec.reserve(data.size());

	for (auto element : data)
	{
		particle_vec.push_back(element.first);
		pixel_vec.push_back(element.second);
	}

	return { particle_vec, pixel_vec };
}

#endif // __LOAD_PRI_FILE_H_

#ifndef __LOAD_PRI_FILE_H_
#define __LOAD_PRI_FILE_H_

#include <vector>
#include <tuple>
#include <fstream>
#include <string>
#include "../core/particle.h"

#if !CUDA_AVAILABLE
struct int2 { int x, y; };
#endif

std::vector<std::pair<particle, int2>> load_pri_file(std::string const & filename)
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
		ifs.read(reinterpret_cast<char*>(&pixel), sizeof(pixel));

		particle primary
		{
			{ buffer[0], buffer[1], buffer[2] }, // pos
			{ buffer[3], buffer[4], buffer[5] }, // dir
			buffer[6]                            // energy
		};

		if (ifs.eof())
			break;

		particle_vec.push_back({ primary, pixel });

		//if (particle_vec.size() > 100000)
		//	break;
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
std::pair<std::vector<particle>, std::vector<int2>> separate_pairs(std::vector<std::pair<particle, int2>> data)
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

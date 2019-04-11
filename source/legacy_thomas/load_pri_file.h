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

size_t num_pris(std::string const & filename)
{
	// Determine number of electrons.
	// We don't use C++17 std::filesystem because we want to be C++11 compatible.
	// Officially, tellg() is not guaranteed to give us the file size, but that
	// is OK because in practice it always does that and it's only indicative.
	std::ifstream ifs(filename, std::ifstream::binary | std::ifstream::ate);
	if (!ifs.is_open())
		return 0;

	size_t size = ifs.tellg();
	return size / (7*sizeof(float) + 2*sizeof(int));
}

std::pair<std::vector<particle>, std::vector<int2>> load_pri_file(std::string const & filename, vec3 min_pos, vec3 max_pos)
{
	std::vector<particle> particle_vec;
	std::vector<int2> pixel_vec;

	// Keep track of the number of times a primary electron is wrong,
	// and keep an example.
	std::pair<size_t, vec3> outside_geom;
	std::pair<size_t, float> low_energy;
	std::pair<size_t, float> high_energy;
	std::pair<size_t, vec3> direction;

	// Reserve expected amount of memory required.
	const size_t estimated_number = num_pris(filename);
	particle_vec.reserve(estimated_number);
	pixel_vec.reserve(estimated_number);

	std::ifstream ifs(filename, std::ifstream::binary);
	if (!ifs.is_open())
		return { particle_vec, pixel_vec };

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
		{
			if (outside_geom.first == 0)
				outside_geom.second = primary.pos;
			++outside_geom.first;
		}

		if (!(primary.kin_energy > EPSILON))
		{
			if (low_energy.first == 0)
				low_energy.second = primary.kin_energy;
			++low_energy.first;
		}

		if (primary.kin_energy > K_max)
		{
			if (high_energy.first == 0)
				high_energy.second = primary.kin_energy;
			++high_energy.first;
		}

		if (!(std::abs(primary.dir.x) > EPSILON || std::abs(primary.dir.y) > EPSILON || std::abs(primary.dir.z) > EPSILON))
		{
			if (direction.first == 0)
				direction.second = primary.dir;
			++direction.first;
		}

		particle_vec.push_back(primary);
		pixel_vec.push_back(pixel);
	}
	ifs.close();

	if (outside_geom.first > 0)
		std::clog << "WARNING: " << outside_geom.first
			<< " primary electrons starting outside geometry. Example: ("
			<< outside_geom.second.x << ", "
			<< outside_geom.second.y << ", "
			<< outside_geom.second.z << ") nm." << std::endl;
	if (low_energy.first > 0)
		std::clog << "WARNING: " << low_energy.first
			<< " primary electrons starting with low kinetic energy. Example: "
			<< low_energy.second << " eV." << std::endl;
	if (high_energy.first > 0)
		std::clog << "WARNING: " << high_energy.first
			<< " primary electrons starting with kinetic energy higher than the hard-coded limit."
			<< " Limit is " << (K_max/1000) << " keV, example energy is "
			<< (high_energy.second/1000) << " keV. Complain to the developers." << std::endl;
	if (direction.first > 0)
		std::clog << "WARNING: " << direction.first
			<< " primary electrons start with unphysical direction vector. Example: ("
			<< direction.second.x << ", "
			<< direction.second.y << ", "
			<< direction.second.z << ")." << std::endl;


	return { particle_vec, pixel_vec };
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


// Helper functions for sorting two vectors simultaneously.
// It does this by finding a sort permutation, which is by no means elegant.
// Given that we don't usually want to sort, we don't care.
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    std::vector<T> const & vec,
    Compare const & compare)
{
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
	return p;
}
template <typename T>
std::vector<T> apply_permutation(
    std::vector<T> const & vec,
    std::vector<size_t> const & p)
{
	std::vector<T> sorted_vec(vec.size());
	std::transform(p.begin(), p.end(), sorted_vec.begin(),
		[&](std::size_t i){ return vec[i]; });
	return sorted_vec;
}

/*
 * Shuffle the primary particles such that the first `prescan_size` are uniformly
 * sampled. The others are left mostly untouched.
 */
void prescan_shuffle(
	std::vector<particle>& particle_vec,
	std::vector<int2>& pixel_vec,
	size_t prescan_size)
{
	std::default_random_engine rng;
	for (size_t i1 = 0; i1 < std::min(prescan_size, particle_vec.size()); ++i1)
	{
		std::uniform_int_distribution<size_t> dist(i1, particle_vec.size()-1);
		const auto i2 = dist(rng);
		std::swap(particle_vec[i1], particle_vec[i2]);
		std::swap(pixel_vec[i1], pixel_vec[i2]);
	}
}

/*
 * The first prescan_size particles are random, uniformly shuffled from the exposure
 * for a well-sampled prescan.
 * The other ones are sorted by Morton index, to make sure that nearby particles in
 * the exposure are also spatially close. This gives a significant speedup in the
 * GPU collision detection routine, but in practice, the gains are usually less
 * than the cost of sorting, especially if the input is already partially sorted.
 */
void sort_pri_file(
	std::vector<particle>& particle_vec,
	std::vector<int2>& pixel_vec)
{
	// Find min, max coordinates of all primaries
	vec3 vmin = particle_vec[0].pos;
	vec3 vmax = vmin;
	for (const auto& p : particle_vec)
	{
		vmin = {
			std::min(vmin.x, p.pos.x),
			std::min(vmin.y, p.pos.y),
			std::min(vmin.z, p.pos.z),
		};
		vmax = {
			std::max(vmax.x, p.pos.x),
			std::max(vmax.y, p.pos.y),
			std::max(vmax.z, p.pos.z),
		};
	}

	// Sort by Morton ordering
	const vec3 vsize = vmax - vmin;
	auto permutation = sort_permutation(particle_vec,
		[vmin, vsize](particle const & p1, particle const & p2) -> bool
	{
		const vec3 p1n = (p1.pos - vmin);
		const vec3 p2n = (p2.pos - vmin);

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

	particle_vec = apply_permutation(particle_vec, permutation);
	pixel_vec = apply_permutation(pixel_vec, permutation);
}

#endif // __LOAD_PRI_FILE_H_

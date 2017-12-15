#ifndef __INDEX_HELPER_H_
#define __INDEX_HELPER_H_

/*
 * Class to assist in computing linear indices for N-dimensional arrays.
 * Indexing is in C order. If the table dimensions are (d1, d2, .., dN), then
 * index (i1, i2, .., iN) maps to (i1*d2*..*dN + i2*d3*..*dN + .. + iN) in
 * linear space.
 */

#include <array>
#include "../../common/variadic.h"

namespace nbl { namespace nd_array { namespace detail {

template<size_t N>
struct index_helper
{
	// Construct given the dimensions (d1, .., dN) of the table
	template<typename... dims>
	index_helper(dims... dimensions)
	{
		static_assert(sizeof...(dimensions) == N, "Wrong number of dimensions");
		set_pitch(dimensions...);
	}

	// Get linear index from (i1, .., iN)
	template<typename... idxs>
	size_t operator()(idxs... indices) const
	{
		static_assert(sizeof...(indices) == N, "Wrong number of indices");
		return get_index(indices...);
	}


	// std::array<size_t, N> version of operator() above.
	size_t operator()(std::array<size_t, N> const & indices) const
	{
		return get_index_array(indices, make_index_sequence<N>());
	}

private:
	size_t pitches[N-1];


	template<typename... idxs>
	size_t get_index(size_t first, idxs... rest) const
	{
		constexpr size_t current_element = N - sizeof...(rest) - 1;
		return first * pitches[current_element] + get_index(rest...);
	}
	size_t get_index(size_t first) const
	{
		return first;
	}

	template<typename... dims>
	size_t set_pitch(size_t /*first*/, size_t second, dims... rest)
	{
		constexpr size_t current_element = N - sizeof...(rest) - 2;
		return pitches[current_element] = second * set_pitch(second, rest...);
	}
	size_t set_pitch(size_t /*first*/)
	{
		return 1;
	}


	template<size_t... s>
	size_t get_index_array(std::array<size_t, N> const & indices, index_sequence<s...>) const
	{
		return get_index(std::get<s>(indices)...);
	}
};


/*
 * Specialization for N=1, to avoid having a size 0 member array.
 */
template<>
struct index_helper<1>
{
	// Construct given the dimensions (d1, .., dN) of the table
	index_helper(size_t)
	{}

	// Get linear index from (i1, .., iN)
	size_t operator()(size_t idx) const
	{
		return idx;
	}

	// std::array<size_t, N> version of operator() above.
	size_t operator()(std::array<size_t, 1> const & indices) const
	{
		return indices[0];
	}
};

}}} // namespace nbl::nd_array::detail

#endif // __INDEX_HELPER_H_

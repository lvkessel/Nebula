#ifndef __INDEX_HELPER_H_
#define __INDEX_HELPER_H_

#include <array>
#include "../../common/variadic.h"

namespace nbl { namespace nd_array { namespace detail {

/**
 * \brief Class to assist in computing linear indices for N-dimensional arrays.
 *
 * Indexing is in C order. If the table dimensions are (d1, d2, .., dN), then
 * index (i1, i2, .., iN) maps to (i1*d2*..*dN + i2*d3*..*dN + .. + iN) in
 * linear space.
 *
 * This class provides a constructor, which needs to know the size of the table
 * along each dimension. The `()` operator then calculates a linear index from
 * the given n-dimensional index.
 *
 * \tparam N Number of dimensions
 */
template<size_t N>
struct index_helper
{
	/**
	 * \brief Construct, given the dimensions `(d1, .., dN)` of the table.
	 *
	 * \param dimensions... Size of the table along each dimension. Size of this
	 *                      parameter pack must be equal to \p N.
	 */
	template<typename... dims>
	index_helper(dims... dimensions)
	{
		static_assert(sizeof...(dimensions) == N, "Wrong number of dimensions");
		set_pitch(dimensions...);
	}

	/**
	 * \brief Get linear index from `N`-dimensional index `(i1, .., iN)`.
	 *
	 * \param indices... "Subscripted" indices. Size of this parameter pack
	 *                   must be equal to \p N.
	 */
	template<typename... idxs>
	size_t operator()(idxs... indices) const
	{
		static_assert(sizeof...(indices) == N, "Wrong number of indices");
		return get_index(indices...);
	}


	/**
	 * \brief Get linear index from `N`-dimensional index `(i1, .., iN)`.
	 *
	 * Version of `operator()` taking a `std::array` instead of `N` separate
	 * parameters.
	 *
	 * \param indices Array of "subscripted" indices.
	 */
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


/**
 * \brief Specialization for `N=1`.
 *
 * Avoids having a size 0 member array in `index_helper<N>`.
 *
 * This class has the same usage as `index_helper<N>`, but its implementation
 * is trivial because a "subscripted" index is the same as a linear index in 1
 * dimension.
 */
template<>
struct index_helper<1>
{
	index_helper(size_t)
	{}

	size_t operator()(size_t idx) const
	{
		return idx;
	}

	size_t operator()(std::array<size_t, 1> const & indices) const
	{
		return indices[0];
	}
};

}}} // namespace nbl::nd_array::detail

#endif // __INDEX_HELPER_H_

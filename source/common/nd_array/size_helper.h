#ifndef __SIZE_HELPER_H_
#define __SIZE_HELPER_H_

#include <numeric>

namespace nbl { namespace nd_array { namespace detail
{

/**
 * \brief Multiplies a list of integers.
 *
 * Used for calculating the total size of an n-dimensional table.
 */
template<typename... dims>
inline size_t size_helper(size_t first, dims... rest)
{
	return first * size_helper(rest...);
}
template<>
inline size_t size_helper(size_t first)
{
	return first;
}


/**
 * \brief Multiplies an array of integers.
 *
 * Used for calculating the total size of an n-dimensional table.
 */
template<size_t N>
inline size_t size_helper(std::array<size_t, N> const & dimensions)
{
	return std::accumulate(dimensions.begin(), dimensions.end(), 1,
		std::multiplies<size_t>());
}

}}} // namespace nbl::nd_array::detail

#endif // __SIZE_HELPER_H_

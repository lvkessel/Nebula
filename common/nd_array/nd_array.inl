#include "size_helper.h"

namespace nbl { namespace nd_array {

template<typename T, size_t N>
template<typename... dims>
nd_array<T, N>::nd_array(dims... dimensions) :
	base_vector_t(detail::size_helper(dimensions...)),
	_indexer(dimensions...),
	_dimensions{dimensions...}
{
	static_assert(sizeof...(dims) == N, "Wrong number of dimensions");
}


template<typename T, size_t N>
size_t nd_array<T, N>::dim(size_t d) const
{
	return _dimensions[d];
}

template<typename T, size_t N>
template<typename... idxs>
T const & nd_array<T, N>::operator()(idxs... indices) const
{
	static_assert(sizeof...(idxs) == N, "Wrong number of indices");
	return base_vector_t::operator[](_indexer(indices...));
}
template<typename T, size_t N>
template<typename... idxs>
T& nd_array<T, N>::operator()(idxs... indices)
{
	static_assert(sizeof...(idxs) == N, "Wrong number of indices");
	return base_vector_t::operator[](_indexer(indices...));
}


template<typename T, size_t N>
template<typename... real_idx_ts>
T nd_array<T, N>::at_linear(real_idx_ts... real_indices) const
{
	static_assert(sizeof...(real_idx_ts) == N, "Wrong number of indices");
	std::array<size_t, N> tmp_arr;
	return at_lin_help(tmp_arr, real_indices...);
}

template<typename T, size_t N>
template<typename... real_idx_ts>
T nd_array<T, N>::at_rounddown(real_idx_ts... real_indices) const
{
	static_assert(sizeof...(real_idx_ts) == N, "Wrong number of indices");
	std::array<size_t, N> tmp_arr;
	return at_rdn_help(tmp_arr, real_indices...);
}


template<typename T, size_t N>
template<typename real_idx_1, typename... real_idx_ts>
T nd_array<T, N>::at_lin_help(std::array<size_t, N>& indices,
	real_idx_1 real_1, real_idx_ts... real_rest) const
{
	constexpr size_t this_dimension = N - sizeof...(real_rest) - 1;
	const size_t low_idx = clamp_index(real_1, 0, dim(this_dimension) - 2);
	const auto frac_idx = real_1 - low_idx;

	// Get data on either side
	indices[this_dimension] = low_idx;
	const T v0 = at_lin_help(indices, real_rest...);
	indices[this_dimension] = low_idx + 1;
	const T v1 = at_lin_help(indices, real_rest...);

	// Return interpolated value
	return (1 - frac_idx) * v0 + frac_idx * v1;
}
template<typename T, size_t N>
T nd_array<T, N>::at_lin_help(std::array<size_t, N>& indices) const
{
	return base_vector_t::operator[](_indexer(indices));
}

template<typename T, size_t N>
template<typename real_idx_1, typename... real_idx_ts>
T nd_array<T, N>::at_rdn_help(std::array<size_t, N>& indices,
	real_idx_1 real_1, real_idx_ts... real_rest) const
{
	constexpr size_t this_dimension = N - sizeof...(real_rest) - 1;
	const size_t low_idx = clamp_index(real_1, 0, dim(this_dimension) - 1);

	indices[this_dimension] = low_idx;
	return at_rdn_help(indices, real_rest...);
}
template<typename T, size_t N>
T nd_array<T, N>::at_rdn_help(std::array<size_t, N>& indices) const
{
	return base_vector_t::operator[](_indexer(indices));
}

// Round real_index down to integer, clamp to range [low, high]
template<typename T, size_t N>
template<typename real_t>
size_t nd_array<T, N>::clamp_index(real_t real_index, size_t low, size_t high)
{
	return static_cast<size_t>(
		std::max(real_t(low), std::min(
			real_index, real_t(high))));
}

}} // namespace nbl::nd_array

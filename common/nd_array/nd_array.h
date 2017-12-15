#ifndef __ND_ARRAY_H_
#define __ND_ARRAY_H_

/*
 * Generic N-dimensional array.
 * Data can be accessed directly with operator(i,j, ...). The functions
 * at_linear(x,y, ...) and at_rounddown(x,y, ...) provide interpolated values
 * at fractional indices. Data is stored in a std::vector under the hood,
 * direct access to this storage is possible using data() or begin()/end().
 * Data is stored contiguously, following the C indexing convention (the last
 * index is the fastest changing one).
 */

#include <array>
#include <vector>
#include "index_helper.h"

namespace nbl { namespace nd_array {

template<typename T, size_t N>
class nd_array
	: private std::vector<T>
{
public:
	using base_vector_t = std::vector<T>;

	// Variadic constructor, takes length along each dimension
	template<typename... dims>
	nd_array(dims... dimensions);

	// Get length along dimension d
	size_t dim(size_t d) const;
	// Get total capacity (product of dim(i))
	using base_vector_t::size;

	// Direct access to an element, unchecked bounds
	template<typename... idxs>
	T const & operator()(idxs... indices) const;
	template<typename... idxs>
	T& operator()(idxs... indices);

	// Linear iterators
	using base_vector_t::begin;
	using base_vector_t::end;


	// Get at a certain index, performing multi-linear interpolation if the
	// indices are floating-point numbers. Extrapolates if out of range.
	template<typename... real_idx_ts>
	T at_linear(real_idx_ts... real_indices) const;
	// Similar, but rounding down to the largest element stored below each
	// index. Rounds up if index is below range.
	template<typename... real_idx_ts>
	T at_rounddown(real_idx_ts... real_indices) const;

	// Direct access to data.
	// Data is in C order: (i, j, k) is at data[i*dim(1)*dim(2) + j*dim(2) + k]
	using base_vector_t::data;

protected:
	// Round real_index down to integer, clamp to range [low, high]
	template<typename real_t>
	static size_t clamp_index(real_t real_index, size_t low, size_t high);

private:
	detail::index_helper<N> _indexer;
	std::array<size_t, N> _dimensions;

	// For use in at_linear and at_rounddown
	template<typename real_idx_1, typename... real_idx_ts>
	T at_lin_help(std::array<size_t, N>&, real_idx_1, real_idx_ts...) const;
	T at_lin_help(std::array<size_t, N>&) const;

	template<typename real_idx_1, typename... real_idx_ts>
	T at_rdn_help(std::array<size_t, N>&, real_idx_1, real_idx_ts...) const;
	T at_rdn_help(std::array<size_t, N>&) const;
};

}} // namespace nbl::nd_array

#include "nd_array.inl"

#endif // __ND_ARRAY_H_

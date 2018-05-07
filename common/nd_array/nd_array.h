#ifndef __ND_ARRAY_H_
#define __ND_ARRAY_H_

#include <array>
#include <vector>
#include "index_helper.h"

namespace nbl { namespace nd_array {

/**
 * \brief Generic N-dimensional array.
 *
 * Data can be accessed directly with `operator(i,j, ...)`. The functions
 * `at_linear(x,y, ...)` and at_rounddown(x,y, ...) provide interpolated values
 * at fractional indices.
 *
 * Data is stored in a std::vector under the hood. Direct access to this storage
 * is possible using `data()` or `begin()`/`end()`. Data is stored contiguously,
 * following the C indexing convention (the last index is the fastest changing
 * one).
 *
 * \tparam T Data type used for storage
 * \tparam N Number of dimensions
 */
template<typename T, size_t N>
class nd_array
	: private std::vector<T>
{
public:
	using base_vector_t = std::vector<T>;

	/**
	 * \brief Variadic constructor, takes length along each dimension.
	 *
	 * \param dimensions... Length along each dimension. Size of this parameter
	 *                      pack must be equal to `N`.
	 */
	template<typename... dims>
	nd_array(dims... dimensions);

	/**
	 * \brief Get size along dimension d
	 */
	size_t dim(size_t d) const;
	/**
	 * \brief Get total capacity (product of dim(i))
	 */
	using base_vector_t::size;

	/**
	 * \brief Direct read-only access to an element.
	 *
	 * It is not verified that the requested data point is actually in range!
	 *
	 * \param indices... Indices to get the element. Size of this parameter pack
	 *                   must be equal to `N`.
	 */
	template<typename... idxs>
	T const & operator()(idxs... indices) const;
	/**
	 * \brief Direct read-write access to an element.
	 *
	 * It is not verified that the requested data point is actually in range!
	 *
	 * \param indices... Indices to get the element at. Size of this parameter
	 *                   pack must be equal to `N`.
	 */
	template<typename... idxs>
	T& operator()(idxs... indices);

	/**
	 * \brief Linear iterator to the beginning of the data.
	 */
	using base_vector_t::begin;
	/**
	 * \brief Linear iterator to one-past-end of the data.
	 */
	using base_vector_t::end;


	/**
	 * \brief Get value with linear interpolation.
	 *
	 * Get read-only acces at a certain index, performing multi-linear
	 * interpolation if the indices are floating-point numbers. Extrapolates if
	 * out of range.
	 *
	 * \param real_indices... Indices to get the element at, may be
	 *                        floating-point numbers. Size of this parameter
	 *                        pack must be equal to `N`.
	 */
	template<typename... real_idx_ts>
	T at_linear(real_idx_ts... real_indices) const;
	/**
	 * \brief Get value, rounding down floating-point indices.
	 *
	 * Get read-only acces at a certain index, rounding indices down to the
	 * nearest integer if they are floating-point numbers. Extrapolates if out
	 * of range.
	 *
	 * \param real_indices... Indices to get the element at, may be
	 *                        floating-point numbers. Size of this parameter
	 *                        pack must be equal to `N`.
	 */
	template<typename... real_idx_ts>
	T at_rounddown(real_idx_ts... real_indices) const;

	/**
	 * \brief Direct access to data.
	 *
	 * Data is in C order. For example, index `(i, j, k)` is at
	 * `data[i*dim(1)*dim(2) + j*dim(2) + k]`.
	 */
	using base_vector_t::data;

protected:
	/**
	 * \brief Round an index down to an integer, clamping to a range.
	 *
	 * \param real_index Index to be rounded down
	 * \param low        Lower boundary of the range to clamp to.
	 * \param high       Upper boundary of the range to clamp to.
	 */
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

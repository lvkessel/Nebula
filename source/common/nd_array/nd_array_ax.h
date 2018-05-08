#ifndef __ND_ARRAY_AX_H_
#define __ND_ARRAY_AX_H_

#include <tuple>
#include "../../common/variadic.h"
#include "nd_array.h"

namespace nbl { namespace nd_array {

/**
 * \brief Simple ::nbl::nd_array::nd_array with additional unit information,
 * plus a scale along each axis.
 *
 * Data from a regular nd_array can only be accessed if the index is known. This
 * class abstracts that and lets the user access data according to a
 * physically meaningful quantity. This is called the "axis" along each
 * dimension.
 *
 * In addition to providing meaningful indices, this table also facilitates in
 * giving the stored data a physically meaningful unit.
 *
 * Specifically for 1D data, this class provides log-log interpolation as well
 * as the usual multi-linear interpolation.
 *
 * \tparam T          Datatype of the data stored
 * \tparam axis_ts... Class representing the scale along each axis. Currently,
 *                    ::nbl::nd_array::ax_list is the only option. Size of this
 *                    parameter pack is equal to the number of dimensions of the
 *                    table.
 */
template<typename T, typename... axis_ts>
class nd_array_ax
	: public nd_array<T, sizeof...(axis_ts)>
{
public:
	static constexpr size_t N = sizeof...(axis_ts); ///< Number of dimensions
	using base_t = nd_array<T, N>;
	using data_type = T;                            ///< Data type for storage
	using value_type = units::quantity<T>;          ///< Data type for the units

	/**
	 * \brief Constructor.
	 *
	 * The size of the underlying `nd_array` is taken from the size of each axis.
	 *
	 * \param axes... The axes to be used.
	 * \param units   Units associated to the data to be stored.
	 */
	nd_array_ax(axis_ts... axes, value_type units)
		: base_t(axes.size()...), _axes(axes...), unit(units)
	{}
	/**
	 * \brief Constructor, using dimensionless units by default.
	 *
	 * TODO: One would normally provide units as a default parameter to the
	 * other constructor and omit this definition. Doesn't seem to work with
	 * with common/variadic.h's ::nbl::make_from_tuple though.
	 */
	nd_array_ax(axis_ts... axes)
		: base_t(axes.size()...), _axes(axes...), unit(units::dimensionless)
	{}

	/**
	 * \brief Get read access to an axis.
	 *
	 * \tparam dim Dimension to get the axis for.
	 */
	template<size_t dim>
	type_at_index<dim, axis_ts...> const & get_axis() const
	{
		return std::get<dim>(_axes);
	}

	/**
	 * \brief Get value with multi-linear interpolation.
	 *
	 * \param coordinates Coordinates (in axis units) to get the value at.
	 */
	template<typename... x_ts>
	value_type get_linear(x_ts... coordinates) const
	{
		return get_lin_helper(make_index_sequence<N>(), coordinates...) * unit;
	}

	/**
	 * \brief Get value, rounding down to nearest stored index.
	 *
	 * \param coordinates Coordinates (in axis units) to get the value at.
	 */
	template<typename... x_ts>
	value_type get_rounddown(x_ts... coordinates) const
	{
		return get_rdn_helper(make_index_sequence<N>(), coordinates...) * unit;
	}

	/**
	 * \brief Log-log interpolation, specifically if there is only one axis.
	 *
	 * \param x Coordinate, in axis units.
	 */
	template<typename x_t>
	value_type get_loglog(x_t x) const
	{
		static_assert(N == 1,
			"log-log interpolation only implemented for 1D arrays.");
		return get_loglog_helper(x) * unit;
	}

	/**
	 * \brief Unit of the stored data.
	 *
	 * The actual data are stored as numbers, without physical units. Their unit
	 * is stored separately in this variable, which acts as a multiplier when
	 * using the `get_xyz()` functions.
	 */
	value_type unit;

private:
	std::tuple<axis_ts...> _axes;

	template<typename... x_ts, size_t... s>
	T get_lin_helper(index_sequence<s...>, x_ts... coordinates) const
	{
		return nd_array<T, N>::at_linear(std::get<s>(_axes).find(coordinates)...);
	}
	template<typename... x_ts, size_t... s>
	T get_rdn_helper(index_sequence<s...>, x_ts... coordinates) const
	{
		return nd_array<T, N>::at_rounddown(std::get<s>(_axes).find(coordinates)...);
	}
	template<typename x_t>
	T get_loglog_helper(x_t coordinate) const
	{
		auto const & x_axis = get_axis<0>();

		const T true_idx = x_axis.find(coordinate);
		const size_t low_idx = base_t::clamp_index(true_idx, 0, base_t::dim(0) - 2);

		const T frac_idx = std::log(coordinate / x_axis[low_idx]) /
			std::log(x_axis[low_idx + 1] / x_axis[low_idx]);
		const T low_value = std::log((*this)(low_idx));
		const T high_value = std::log((*this)(low_idx + 1));

		return std::exp((1 - frac_idx)*low_value + frac_idx*high_value);
	}
};

}} // namespace nbl::nd_array

#endif

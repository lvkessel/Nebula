#ifndef __ND_ARRAY_AX_H_
#define __ND_ARRAY_AX_H_

/*
 * A simple nd_array with additional unit information, plus a scale along each
 * axis. Very useful if the user does not want to know about the discrete
 * nature of the underlying data.
 * 
 * Specifically for 1D data, this class provides log-log interpolation as well
 * as the usual multi-linear interpolation.
 */

#include <tuple>
#include "../../common/variadic.h"
#include "nd_array.h"

namespace nbl { namespace nd_array {

template<typename T, typename... axis_ts>
class nd_array_ax
	: public nd_array<T, sizeof...(axis_ts)>
{
public:
	static constexpr size_t N = sizeof...(axis_ts); // Number of dimensions
	using base_t = nd_array<T, N>;
	using data_type = T;
	using value_type = units::quantity<T>;

	// Construct, given axes.
	// The dimensions of the underlying nd_array is taken from the size of each axis.
	nd_array_ax(axis_ts... axes, value_type units)
		: base_t(axes.size()...), _axes(axes...), _units(units)
	{}
	// TODO: one would normally provide units as a default parameter and omit the definition below.
	// However, Visual Studio does not seem to like that combination with common/variadic.h's construct_unroll.
	nd_array_ax(axis_ts... axes)
		: base_t(axes.size()...), _axes(axes...), _units(units::dimensionless)
	{}

	// Get read access to an axis
	template<size_t dim>
	type_at_index<dim, axis_ts...> const & get_axis() const
	{
		return std::get<dim>(_axes);
	}

	// Perform multi-linear interpolation
	template<typename... x_ts>
	value_type get_linear(x_ts... coordinates) const
	{
		return get_lin_helper(make_index_sequence<N>(), coordinates...) * _units;
	}
	// Same, but rounding down to the largest element below each index.
	template<typename... x_ts>
	value_type get_rounddown(x_ts... coordinates) const
	{
		return get_rdn_helper(make_index_sequence<N>(), coordinates...) * _units;
	}

	// Log-log interpolation, specifically if there is only one axis
	template<typename x_t>
	value_type get_loglog(x_t x) const
	{
		static_assert(sizeof...(axis_ts) == 1, "log-log interpolation only implemented for 1D arrays.");
		return get_loglog_helper(x) * _units;
	}

	// Set the units. This does not convert the underlying data! It merely
	// assigns a different meaning to the numerical values.
	void set_units(value_type units)
	{
		_units = units;
	}

private:
	std::tuple<axis_ts...> _axes;
	units::quantity<data_type> _units;

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

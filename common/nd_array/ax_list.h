#ifndef __AX_LIST_H_
#define __AX_LIST_H_

/*
 * Axis representation as array of consecutive values.
 * Data points are assumed to be consecutive, this is never explicitly checked.
 */

#include <vector>
#include <algorithm>
#include "nd_array.h"
#include "../units/quantity.h"

namespace nbl { namespace nd_array {

template<typename T>
class ax_list : private std::vector<T>
{
public:
	using data_type = T;
	using value_type = units::quantity<T>;
	using base_type = std::vector<T>;

	using base_type::size;

	ax_list() = default;
	ax_list(std::vector<T> const & data,
		value_type units = units::dimensionless)
		: base_type(data), _units(units)
	{}
	ax_list(nd_array<T, 1> const & data,
		value_type units = units::dimensionless)
		: base_type(data.begin(), data.end()), _units(units)
	{}
	ax_list(size_t size,
		value_type units = units::dimensionless)
		: base_type(size), _units(units)
	{}

	value_type operator[](size_t i) const
	{
		return base_type::operator[](i) * _units;
	}

	// Find the position of x along this axis.
	// Returns fractional index if x is not found directly.
	data_type find(value_type x) const
	{
		const data_type xx = x / _units;

		// Estimate true index, even if out of range.
		const auto high_iterator = std::lower_bound(base_type::begin(),
			base_type::end(), xx);

		// Not using std::clamp, so we won't require more than C++11
		const size_t high_index = std::max<size_t>(1, std::min(
			std::distance(base_type::begin(), high_iterator),
			static_cast<typename base_type::difference_type>(size()) - 1));
		const data_type high_value = base_type::operator[](high_index);
		const data_type low_value = base_type::operator[](high_index - 1);

		return high_index +
			(xx - high_value) / (high_value - low_value);
	}

	value_type units() const
	{
		return _units;
	}

private:
	units::quantity<data_type> _units;
};

}} // namespace nbl::nd_array

#endif

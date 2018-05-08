#ifndef __AX_LIST_H_
#define __AX_LIST_H_

#include <vector>
#include <algorithm>
#include "nd_array.h"
#include "../units/quantity.h"

namespace nbl { namespace nd_array {

/**
 * \brief Axis representation as array of consecutive values.
 *
 * Data points are assumed to be ordered strictly ascending, but this is never
 * explicitly checked.
 *
 * This class carries units (represented by ::nbl::units::quantity) as well as
 * values. Each data point is stored as a scalar, the unit (the same for all
 * data points) is stored separately. The scalar data point is then multiplied
 * by the unit when data from this class is requested.
 *
 * \tparam T Data type used to store the points.
 */
template<typename T>
class ax_list : private std::vector<T>
{
public:
	using data_type = T;                   ///< Data type for storage
	using value_type = units::quantity<T>; ///< Data type for the units.
	using base_type = std::vector<T>;

	using base_type::size;

	/// Default constructor. Initializes to no values, unit is undefined.
	ax_list() = default;

	/**
	 * \brief Construct with initial data points and units.
	 *
	 * \param data  Initial data points.
	 * \param units Units associated with this axis.
	 */
	ax_list(std::vector<T> const & data,
		value_type units = units::dimensionless)
		: base_type(data), _units(units)
	{}
	/**
	 * \brief Construct with initial data points and units.
	 *
	 * \param data  Initial data points.
	 * \param units Units associated with this axis.
	 */
	ax_list(nd_array<T, 1> const & data,
		value_type units = units::dimensionless)
		: base_type(data.begin(), data.end()), _units(units)
	{}
	/**
	 * \brief Construct with initial size and units, but without data.
	 *
	 * Data points are zero-initialized.
	 *
	 * \param size  Initial data size.
	 * \param units Units associated with this axis.
	 */
	ax_list(size_t size,
		value_type units = units::dimensionless)
		: base_type(size), _units(units)
	{}

	/**
	 * \brief Access a value stored at a given index.
	 *
	 * Note: this is read-only access. This function multiplies the data point
	 * with the unit.
	 *
	 * \param i Index to be accessed.
	 */
	value_type operator[](size_t i) const
	{
		return base_type::operator[](i) * _units;
	}

	/**
	 * \brief Find the position of a certain value along this axis.
	 *
	 * Returns fractional index if the value is not found directly.
	 *
	 * \param x Value to be found, including the correct units.
	 */
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

	/// Read access to the units.
	value_type units() const
	{
		return _units;
	}

private:
	units::quantity<data_type> _units;
};

}} // namespace nbl::nd_array

#endif

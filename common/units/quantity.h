#ifndef __QUANTITY_H_
#define __QUANTITY_H_

/*
 * A quantity is a number plus a dimension.
 */

#include <cmath>
#include <stdexcept>
#include "dimension.h"

namespace nbl { namespace units {

template<typename T>
struct quantity
{
	using value_type = T;

	value_type value;
	dimension units;

	bool dimensionless() const
	{
		return units == dimensions::dimensionless;
	}

	operator T() const
	{
		if (!dimensionless())
			throw std::runtime_error("Can only convert dimensionless quantity"
				" to scalar");
		return value;
	}

	quantity& operator+=(quantity const & rhs)
	{
		if (units != rhs.units)
			throw std::runtime_error("Adding incompatible units.");
		value += rhs.value;
		return *this;
	}
	quantity& operator-=(quantity const & rhs)
	{
		if (units != rhs.units)
			throw std::runtime_error("Adding incompatible units.");
		value -= rhs.value;
		return *this;
	}
	quantity& operator*=(quantity const & rhs)
	{
		value *= rhs.value;
		units *= rhs.units;
		return *this;
	}
	quantity& operator/=(quantity const & rhs)
	{
		value /= rhs.value;
		units /= rhs.units;
		return *this;
	}

	quantity& operator*=(value_type const rhs)
	{
		value *= rhs;
		return *this;
	}
	quantity& operator/=(value_type const rhs)
	{
		value /= rhs;
		return *this;
	}

	const quantity operator+(quantity const & rhs) const
	{
		return quantity(*this) += rhs;
	}
	const quantity operator-(quantity const & rhs) const
	{
		return quantity(*this) -= rhs;
	}
	const quantity operator*(quantity const & rhs) const
	{
		return quantity(*this) *= rhs;
	}
	const quantity operator/(quantity const & rhs) const
	{
		return quantity(*this) /= rhs;
	}
	const quantity operator*(value_type const rhs) const
	{
		return quantity(*this) *= rhs;
	}
	const quantity operator/(value_type const rhs) const
	{
		return quantity(*this) /= rhs;
	}

	bool operator==(quantity const & rhs) const
	{
		return value == rhs.value
			&& units == rhs.units;
	}
	bool operator!=(quantity const & rhs) const
	{
		return !(*this == rhs);
	}
};

template<typename T>
inline quantity<T> pow(quantity<T> q, int power)
{
	return
	{
		std::pow(q.value, power),
		pow(q.units, power)
	};
}

template<typename T>
inline const quantity<T> operator*(T const v, quantity<T> const & q)
{
	return q * v;
}
template<typename T>
inline const quantity<T> operator/(T const v, quantity<T> const & q)
{
	return
	{
		v / q.value,
		dimensions::dimensionless / q.units
	};
}

}} // namespace nbl::units

#endif

#ifndef __HDF5_FILE_H_
#define __HDF5_FILE_H_

/*
 * Simple interface to HDF5 files as we use them for storing material data.
 * We can read:
 *   - Properties
 *     Stored using the HDF5 attribute system, these can either be text or a
 *     number + unit.
 *   - N-dimensional data sets, with axis data stored using the HDF5 dimension
 *     scale API. The data set and axes may have a "unit" attribute. Such data
 *     is returned as an nd_array_ax<double, ...> with N ax_lists as axes.
 */

#include <string>
#include <H5Cpp.h>
#include "../common/units/unit_system.h"
#include "../common/nd_array/nd_array.h"
#include "../common/nd_array/nd_array_ax.h"
#include "../common/nd_array/ax_list.h"

namespace nbl {

class hdf5_file
{
public:
	// Define the hdf5_table<N> data type, which is
	// nd_array_ax<double, ax_list<double>, ax_list<double>, ...>
	template<typename... axes>
	using nd_array_double = nd_array::nd_array_ax<double, axes...>;
	template<size_t N>
	using hdf5_table = repeat<nd_array::ax_list<double>, N, nd_array_double>;


	// Constructor: open a file for reading
	hdf5_file(std::string const & filename);

	// Get property, either a string or a quantity
	std::string get_property_string(std::string const & name) const;
	units::quantity<double> get_property_quantity(std::string const & name,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

	// Get an N-dimensional dataset, with units and axis information.
	// If axis information is not given in the file, it is a dimensionless
	// linear axis from 0 to 1 (inclusive)
	template<size_t N>
	hdf5_table<N> get_table_axes(std::string const & name,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

private:
	H5::H5File _file;
};

} // namespace nbl

#endif // __H5_FILE_H_

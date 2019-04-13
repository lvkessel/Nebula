#ifndef __HDF5_FILE_H_
#define __HDF5_FILE_H_

#include <string>
#include <hdf5.h>
#include "../common/units/unit_system.h"
#include "../common/nd_array/nd_array.h"
#include "../common/nd_array/nd_array_ax.h"
#include "../common/nd_array/ax_list.h"

namespace nbl {

/**
 * \brief Simple interface to HDF5 files as we use them for storing material data.
 *
 * We can read:
 *   - Properties.
 *     Stored using the HDF5 attribute system, these can either be text or a
 *     number + unit.
 *   - N-dimensional data sets, with axis data stored using the HDF5 dimension
 *     scale API. The data set and axes may have a "unit" attribute. Such data
 *     is returned as an ::nbl::nd_array::nd_array_ax <double, ...>, where `...`
 *     is `N` times ::nbl::nd_array::ax_list as axes.
 */
class hdf5_file
{
public:
	template<typename... axes>
	using nd_array_double = nd_array::nd_array_ax<double, axes...>;
	/**
	 * \brief The `N`-dimensional array dataype that can be read from a HDF5 file.
	 *
	 * Explicitly, this is `nd_array_ax<double, ax_list<double>, ax_list<double>, ...>`
	 */
	template<size_t N>
	using hdf5_table = repeat<nd_array::ax_list<double>, N, nd_array_double>;


	/**
	 * \brief Constructor: open a file for reading
	 */
	hdf5_file(std::string const & filename);

	/**
	 * \brief Destructor.
	 */
	~hdf5_file();

	/**
	 * \brief Get a string property. Throws `std::runtime_error` exception if
	 * not found.
	 *
	 * \param name Name of the property to be found.
	 */
	std::string get_property_string(std::string const & name) const;
	/**
	 * \brief Get a quantity property. Throws `std::runtime_error` exception if
	 * not found.
	 *
	 * \param name   Name of the property to be found.
	 * \param parser Class used to parse the quantity's units
	 */
	units::quantity<double> get_property_quantity(std::string const & name,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

	/**
	 * \brief Get a string property. Returns default if not found.
	 *
	 * \param name     Name of the property to be found
	 * \param _default Value returned if property is not found
	 */
	std::string get_property_string(std::string const & name,
		std::string const & _default) const;
	/**
	 * \brief Get a quantity property. Returns default if not found.
	 *
	 * \param name     Name of the property to be found
	 * \param _default Value returned if property is not found
	 * \param parser   Class used to parse the quantity's units
	 */
	units::quantity<double> get_property_quantity(std::string const & name,
		units::quantity<double> const & _default,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

	/**
	 * \brief Get an N-dimensional dataset, with units and axis information.
	 *
	 * If axis information is not given in the file, it is a dimensionless
	 * linear axis from 0 to 1 (inclusive).
	 *
	 * Throws `std::runtime_error` if the table is not found or if its dimension
	 * in the file is not equal to `N`.
	 *
	 * \tparam N      Number of dimensions expected
	 * \param  name   Name of the dataset
	 * \param  parser Class used to parse the quantity's units
	 */
	template<size_t N>
	hdf5_table<N> get_table_axes(std::string const & name,
		units::unit_parser<double> const & parser = units::default_unit_parser()) const;

private:
	hid_t _file;
};

} // namespace nbl

#endif // __H5_FILE_H_

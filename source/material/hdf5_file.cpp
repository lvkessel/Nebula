#include "hdf5_file.h"
#include <hdf5_hl.h>

namespace nbl {

namespace
{
	/*
	 * Get size of dataspace in each dimension.
	 * Verifies that the table has the expected number of dimensions.
	 */
	template<size_t N>
	inline std::array<hsize_t, N> get_size(hid_t dataspace)
	{
		if (H5Sget_simple_extent_ndims(dataspace) != N)
			throw std::runtime_error("Dataspace has unexpected dimension");

		std::array<hsize_t, N> dim;
		H5Sget_simple_extent_dims(dataspace, dim.data(), nullptr);
		return dim;
	}


	/*
	 * Get native string data type for use with the HDF5 library
	 */
	inline hid_t string_type()
	{
		hid_t tid = H5Tcopy(H5T_C_S1);
		H5Tset_size(tid, H5T_VARIABLE);
		return tid;
	}


	/*
	 * Read an attribute, expecting a string
	 */
	inline std::string read_attribute_string(hid_t attribute)
	{
		char* buffer;
		hid_t string_t = string_type();

		H5Aread(attribute, string_t, &buffer);
		std::string result(buffer);

		H5Tclose(string_t);
		H5free_memory(buffer);
		return result;
	}

	/*
	 * Read an attribute, expecting a quantity (value + unit)
	 */
	inline units::quantity<double> read_attribute_quantity(
		hid_t attribute, units::unit_parser<double> const & parser)
	{
		struct data_struct
		{
			double value;
			char* unit;
		};

		// Create datatype
		hid_t string_t = string_type();
		hid_t data_type = H5Tcreate(H5T_COMPOUND, sizeof(data_struct));
		H5Tinsert(data_type, "value", HOFFSET(data_struct, value), H5T_NATIVE_DOUBLE);
		H5Tinsert(data_type, "unit", HOFFSET(data_struct, unit), string_t);

		// Read from file
		data_struct h5_result;
		H5Aread(attribute, data_type, &h5_result);

		// Copy to result struct
		const units::quantity<double> result = h5_result.value *
			parser.parse_unit(h5_result.unit);

		// Clean up and return
		H5free_memory(h5_result.unit);
		H5Tclose(data_type);
		H5Tclose(string_t);
		return result;
	}


	/*
	 * Read entire N-dimensional dataset as doubles.
	 * Verifies that the dataset has the required dimensionality.
	 * Returns this as an nd_array<double, N>.
	 */
	template<size_t N>
	inline nd_array::nd_array<double, N> read_dataset(hid_t dataset)
	{
		hid_t dataspace = H5Dget_space(dataset);
		std::array<hsize_t, N> dimensions = get_size<N>(dataspace);
		H5Sclose(dataspace);

		auto table = make_from_tuple<nd_array::nd_array<double,N>>(dimensions);
		H5Dread(dataset, H5T_NATIVE_DOUBLE,
			H5S_ALL, H5S_ALL, H5P_DEFAULT, table.data());

		return table;
	}


	/*
	 * Read dimension scales attached to a dataset.
	 * If multiple dimension scales are found, the first is returned; if none
	 * are found, this function creates one with points evenly spaced between
	 * zero and one (inclusive).
	 */
	template<size_t N>
	inline std::array<nd_array::ax_list<double>, N> get_scales(
		hid_t dataset,
		std::array<hsize_t, N> const & dimensions,
		units::unit_parser<double> const & parser)
	{
		std::array<nd_array::ax_list<double>, N> scales;

		for (size_t dim = 0; dim < N; ++dim)
		{
			// Find out which dimension scale, if any, is attached.
			// If there are multiple, just get the first one.

			// Assemble data to be sent to the iteration function
			std::pair<
				nd_array::ax_list<double>*,
				units::unit_parser<double> const *>
			itdata { &scales[dim], &parser };

			// Iterate through the data scales attached, if any, and retrieve data
			H5DSiterate_scales(dataset, (unsigned int)dim, nullptr,
				[&](hid_t /*did*/, unsigned /*dim*/, hid_t dsid, void* data) -> herr_t
				{
					// Recover itdata
					auto itdata = *reinterpret_cast<std::pair<
						nd_array::ax_list<double>*,
						units::unit_parser<double> const *>*>(data);

					// Read unit associated to dimension scale
					hid_t at_units = H5Aopen(dsid, "units", H5P_DEFAULT);
					auto unit = itdata.second->parse_unit(
						read_attribute_string(at_units));
					H5Aclose(at_units);

					// Read the dataset
					*(itdata.first) = { read_dataset<1>(dsid), unit };

					// Return 1 to stop the iteration
					return 1;
				}, &itdata);


			// If no dimension scale was attached: fill with 0 to 1 (inclusive)
			if (scales[dim].size() == 0)
			{
				std::vector<double> data(dimensions[dim]);
				for (size_t i = 0; i < data.size(); ++i)
					data[i] = double(i) / (data.size()-1);
				scales[dim] = { data };
			}


			if (scales[dim].size() != dimensions[dim])
				throw std::runtime_error("Dimension scale has unexpected size.");
		}

		return scales;
	}

	/*
	 * Read entire N-dimensional dataset as doubles, including dimension scales.
	 * This function is aware of the units belonging to the dataset and scales.
	 * Each dimension without associated dimension scale in the file is assigned
	 * a linear scale between zero and one.
	 */
	template<size_t N>
	inline hdf5_file::hdf5_table<N> read_axis_dataset(hid_t dataset,
		units::unit_parser<double> const & parser)
	{
		// Get the dimensions and scales
		hid_t dataspace = H5Dget_space(dataset);
		const std::array<hsize_t, N> dimensions = get_size<N>(dataspace);
		H5Sclose(dataspace);
		const auto scales = get_scales<N>(dataset, dimensions, parser);

		// Read the actual data + units
		auto table = make_from_tuple<hdf5_file::hdf5_table<N>>(scales);
		H5Dread(dataset, H5T_NATIVE_DOUBLE,
			H5S_ALL, H5S_ALL, H5P_DEFAULT, table.data());

		hid_t at_units = H5Aopen(dataset, "units", H5P_DEFAULT);
		table.unit = parser.parse_unit(read_attribute_string(at_units));
		H5Aclose(at_units);

		return table;
	}
} // anonymous namespace

hdf5_file::hdf5_file(std::string const & filename)
{
	_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (_file < 0)
		throw std::runtime_error("Unable to open file " + filename);
}
hdf5_file::~hdf5_file()
{
	 H5Fclose(_file);
}

std::string hdf5_file::get_property_string(std::string const & name) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		throw std::runtime_error("Could not read string property '" + name + '\'');

	auto result = read_attribute_string(attribute);
	H5Aclose(attribute);
	return result;
}
units::quantity<double> hdf5_file::get_property_quantity(
	std::string const & name,
	units::unit_parser<double> const & parser) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		throw std::runtime_error("Could not read quantity property '" + name + '\'');

	auto result = read_attribute_quantity(attribute, parser);
	H5Aclose(attribute);
	return result;
}

std::string hdf5_file::get_property_string(
	std::string const & name,
	std::string const & _default) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		return _default;

	auto result = read_attribute_string(attribute);
	H5Aclose(attribute);
	return result;
}
units::quantity<double> hdf5_file::get_property_quantity(
	std::string const & name,
	units::quantity<double> const & _default,
	units::unit_parser<double> const & parser) const
{
	hid_t attribute = H5Aopen(_file, name.c_str(), H5P_DEFAULT);
	if (attribute < 0)
		return _default;

	auto result = read_attribute_quantity(attribute, parser);
	H5Aclose(attribute);
	return result;
}

template<size_t N>
hdf5_file::hdf5_table<N> hdf5_file::get_table_axes(
	std::string const & name,
	units::unit_parser<double> const & parser) const
{
	hid_t dataset = H5Dopen(_file, name.c_str(), H5P_DEFAULT);
	if (dataset < 0)
		throw std::runtime_error("Could not read table '" + name + '\'');

	auto result = read_axis_dataset<N>(dataset, parser);
	H5Dclose(dataset);
	return result;
}


// Explicit instantiations
template hdf5_file::hdf5_table<1> hdf5_file::get_table_axes<1>(std::string const &, units::unit_parser<double> const &) const;
template hdf5_file::hdf5_table<2> hdf5_file::get_table_axes<2>(std::string const &, units::unit_parser<double> const &) const;
template hdf5_file::hdf5_table<3> hdf5_file::get_table_axes<3>(std::string const &, units::unit_parser<double> const &) const;

} // namespace nbl

#include <H5Cpp.h>
#include <hdf5_hl.h>
#include "hdf5_file.h"

namespace nbl {

namespace
{
	/*
	 * Get size of dataspace in each dimension.
	 * Verifies that the table has the expected number of dimensions.
	 */
	template<size_t N>
	inline std::array<hsize_t, N> get_size(H5::DataSpace const & dataspace)
	{
		if (dataspace.getSimpleExtentNdims() != N)
			throw std::runtime_error("Dataspace has unexpected dimension");

		std::array<hsize_t, N> dim;
		dataspace.getSimpleExtentDims(dim.data());
		return dim;
	}


	/*
	 * Get native std::string data type for the HDF5 library
	 */
	inline H5::StrType string_type()
	{
		return H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);
	}

	/*
	 * Free temporary (string) data allocated by the HDF5 library.
	 */
	inline void HDFree(void* pointer)
	{
		free(pointer);
	}


	/*
	 * Read an attribute, expecting a string
	 */
	inline std::string read_attribute_string(H5::Attribute const & attribute)
	{
		std::string result;
		attribute.read(string_type(), result);
		return result;
	}

	/*
	 * Read an attribute, expecting a quantity (value + unit)
	 */
	inline units::quantity<double> read_attribute_quantity(
		H5::Attribute const & attribute, units::unit_parser<double> const & parser)
	{
		struct data_struct
		{
			double value;
			char* unit;
		};

		// Read from file
		H5::CompType data_type(sizeof(data_struct));
		data_type.insertMember("value", HOFFSET(data_struct, value),
			H5::PredType::NATIVE_DOUBLE);
		data_type.insertMember("unit", HOFFSET(data_struct, unit),
			string_type());
		data_struct h5_result;
		attribute.read(data_type, &h5_result);

		// Copy to result struct
		const units::quantity<double> result = h5_result.value *
			parser.parse_unit(h5_result.unit);

		// Clean up and return
		HDFree(h5_result.unit);
		return result;
	}


	/*
	 * Read entire N-dimensional dataset as doubles.
	 * Verifies that the dataset has the required dimensionality.
	 * Returns this as an nd_array<double, N>.
	 */
	template<size_t N>
	inline nd_array::nd_array<double, N> read_dataset(H5::DataSet dataset)
	{
		H5::DataSpace dataspace = dataset.getSpace();
		std::array<hsize_t, N> dimensions = get_size<N>(dataspace);

		auto table = construct_unroll<nd_array::nd_array<double,N>>(dimensions);
		dataset.read(table.data(), H5::PredType::NATIVE_DOUBLE);

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
		H5::DataSet const & dataset,
		std::array<hsize_t, N> const & dimensions,
		units::unit_parser<double> const & parser)
	{
		std::array<nd_array::ax_list<double>, N> scales;

		for (size_t dim = 0; dim < N; ++dim)
		{
			// Find out which dimension scale, if any, is attached.
			// If there are multiple, just get the first one.
			std::pair<bool, H5::DataSet> scale_data(false, {});
			H5DSiterate_scales(dataset.getId(), (unsigned int)dim, nullptr,
				[](hid_t /*did*/, unsigned /*dim*/, hid_t dsid, void* data) -> herr_t
				{
					(*reinterpret_cast<std::pair<bool, H5::DataSet>*>(data)) =
						{ true, H5::DataSet(dsid) };
					return 0;
				}, &scale_data);

			if (scale_data.first == true)
			{
				// Read the dimension scale
				auto unit = parser.parse_unit(read_attribute_string(
					scale_data.second.openAttribute("units")));
				scales[dim] = { read_dataset<1>(scale_data.second), unit };
				if (scales[dim].size() != dimensions[dim])
					throw std::runtime_error("Dimension scale has unexpected size.");
			}
			else
			{
				// No dimension scale attached, fill with 0 to 1 (inclusive)
				std::vector<double> data(dimensions[dim]);
				for (size_t i = 0; i < data.size(); ++i)
					data[i] = double(i) / (data.size()-1);
				scales[dim] = { data };
			}
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
	inline hdf5_file::hdf5_table<N> read_axis_dataset(H5::DataSet dataset,
		units::unit_parser<double> const & parser)
	{
		// Get the dimensions and scales
		const H5::DataSpace dataspace = dataset.getSpace();
		const std::array<hsize_t, N> dimensions = get_size<N>(dataspace);
		const auto scales = get_scales<N>(dataset, dimensions, parser);

		// Read the actual data + units
		auto table = construct_unroll<hdf5_file::hdf5_table<N>>(scales);
		dataset.read(table.data(), H5::PredType::NATIVE_DOUBLE);
		table.unit = parser.parse_unit(read_attribute_string(
			dataset.openAttribute("units")));
		return table;
	}
} // anonymous namespace

// The weird (H5::Exception::dontPrint(), filename) construction exists because
// we want to call H5::Exception::dontPrint() before initializing _file. The
// comma operator performs dontPrint and evaluates to filename.
hdf5_file::hdf5_file(std::string const & filename)
	try : _file((H5::Exception::dontPrint(), filename), H5F_ACC_RDONLY)
{}
catch(H5::Exception const &)
{
	throw std::runtime_error("Could not open HDF5 file " + filename);
}

std::string hdf5_file::get_property_string(std::string const & name) const
{
	try
	{
		return read_attribute_string(_file.openAttribute(name));
	}
	catch (H5::Exception const &)
	{
		throw std::runtime_error("Could not read string property '" + name + '\'');
	}
}
units::quantity<double> hdf5_file::get_property_quantity(
	std::string const & name,
	units::unit_parser<double> const & parser) const
{
	try
	{
		return read_attribute_quantity(
			_file.openAttribute(name),
			parser);
	}
	catch (H5::Exception const &)
	{
		throw std::runtime_error("Could not read quantity property '" + name + '\'');
	}
}

std::string hdf5_file::get_property_string(
	std::string const & name,
	std::string const & _default) const
{
	try
	{
		return read_attribute_string(_file.openAttribute(name));
	}
	catch (H5::Exception const &)
	{
		return _default;
	}
}
units::quantity<double> hdf5_file::get_property_quantity(
	std::string const & name,
	units::quantity<double> const & _default,
	units::unit_parser<double> const & parser) const
{
	try
	{
		return read_attribute_quantity(
			_file.openAttribute(name),
			parser);
	}
	catch (H5::Exception const &)
	{
		return _default;
	}
}

template<size_t N>
hdf5_file::hdf5_table<N> hdf5_file::get_table_axes(
	std::string const & name,
	units::unit_parser<double> const & parser) const
{
	try
	{
		return read_axis_dataset<N>(_file.openDataSet(name), parser);
	}
	catch (H5::Exception const &)
	{
		throw std::runtime_error("Could not read table '" + name + '\'');
	}
}


// Explicit instantiations
template hdf5_file::hdf5_table<1> hdf5_file::get_table_axes<1>(std::string const &, units::unit_parser<double> const &) const;
template hdf5_file::hdf5_table<2> hdf5_file::get_table_axes<2>(std::string const &, units::unit_parser<double> const &) const;
template hdf5_file::hdf5_table<3> hdf5_file::get_table_axes<3>(std::string const &, units::unit_parser<double> const &) const;

} // namespace nbl

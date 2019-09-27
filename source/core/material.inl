template<typename... scatter_types>
CPU material<scatter_list<scatter_types...>>::material(material_legacy_thomas const & mat_legacy)
	: scatter_list<scatter_types...>(scatter_types::create(mat_legacy)...),
	barrier(static_cast<real>(mat_legacy.barrier() / constant::ec))
{
}

template<typename... scatter_types>
CPU material<scatter_list<scatter_types...>>::material(nbl::hdf5_file const & mat)
	: scatter_list<scatter_types...>(scatter_types::create(mat)...),
	barrier(static_cast<real>(mat.get_property_quantity("barrier") / nbl::units::eV))
{
}

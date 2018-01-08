#ifndef __MATERIAL_MANAGER_H_
#define __MATERIAL_MANAGER_H_

namespace nbl { namespace drivers {

template<typename material_type, bool gpu_flag>
struct _material_manager_factory;

template<typename material_type, bool gpu_flag>
class material_manager
{
public:
	using material_t = material_type;
	using material_index_t = int;

	static CPU material_manager<material_t, gpu_flag> create(std::vector<material_t> const & materials);
	static CPU void destroy(material_manager & manager);

	// Direct access to a material, no bounds checking
	inline PHYSICS material_t & operator[](material_index_t i);
	inline PHYSICS material_t const & operator[](material_index_t i) const;

	// Checks if the material is not a "special" vacuum-type material (see index_enum).
	// This function only checks if the provided index >= 0, that is, it does not check
	// if the provided index < _capacity.
	// TODO: invent nicer name...
	inline PHYSICS bool is_physical(material_index_t material_idx) const;

	enum index_enum : material_index_t
	{
		NOP = -128,            // TODO: support this
		TERMINATOR = -127,
		DETECTOR = -126,
		DETECTOR_LT50 = -125,
		DETECTOR_GE50 = -124,
		VACUUM = -123,
		MIRROR = -122
	};

private:
	material_index_t _capacity = 0;
	material_t* _materials = nullptr;

	friend struct _material_manager_factory<material_t, gpu_flag>;
};

}} // namespace nbl::drivers

#include "material_manager.inl"

#endif // __MATERIAL_MANAGER_H_

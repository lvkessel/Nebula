#include <cstring>

namespace nbl { namespace util {

template<typename T, bool gpu_flag>
CPU table_3D<T, gpu_flag> table_3D<T, gpu_flag>::create(
	real x_min, real x_max, size_t width,
	real y_min, real y_max, size_t height,
	real z_min, real z_max, size_t depth,
	T* data)
{
	table_3D<T, gpu_flag> table;
	_table_3D_factory<T, gpu_flag>::allocate(table, width, height, depth);
	if (data != nullptr)
		_table_3D_factory<T, gpu_flag>::set(table, data);

	table._x_min = x_min;
	table._x_step = (width - 1) / (x_max - x_min);
	table._y_min = y_min;
	table._y_step = (height - 1) / (y_max - y_min);
	table._z_min = z_min;
	table._z_step = (depth - 1) / (z_max - z_min);

	return table;
}
template<typename T, bool gpu_flag>
CPU void table_3D<T, gpu_flag>::destroy(table_3D<T, gpu_flag> & table)
{
	_table_3D_factory<T, gpu_flag>::free(table);
}

template<typename T, bool gpu_flag>
template<typename callback_function>
CPU void table_3D<T, gpu_flag>::mem_scope(callback_function callback)
{
	_table_3D_factory<T, gpu_flag>::mem_scope(*this, callback);
}

template<typename T, bool gpu_flag>
PHYSICS T & table_3D<T, gpu_flag>::operator()(size_t i, size_t j, size_t k)
{
	return *(reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(_data) + (j + k*_height) * _pitch) + i);
}

template<typename T, bool gpu_flag>
PHYSICS T const & table_3D<T, gpu_flag>::operator()(size_t i, size_t j, size_t k) const
{
	return *(reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(_data) + (j + k*_height) * _pitch) + i);
}

template<typename T, bool gpu_flag>
PHYSICS T table_3D<T, gpu_flag>::get(real x, real y, real z) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;
	const real z_index = (z - _z_min) * _z_step;

	const int low_x = clampi(floorr(x_index), 0, _width - 2);
	const int low_y = clampi(floorr(y_index), 0, _height - 2);
	const int low_z = clampi(floorr(z_index), 0, _depth - 2);

	const T v000 = (*this)(low_x, low_y, low_z);
	const T v100 = (*this)(low_x + 1, low_y, low_z);
	const T v010 = (*this)(low_x, low_y + 1, low_z);
	const T v110 = (*this)(low_x + 1, low_y + 1, low_z);
	const T v001 = (*this)(low_x, low_y, low_z + 1);
	const T v101 = (*this)(low_x + 1, low_y, low_z + 1);
	const T v011 = (*this)(low_x, low_y + 1, low_z + 1);
	const T v111 = (*this)(low_x + 1, low_y + 1, low_z + 1);

	const real frac_x = x_index - low_x;
	const real frac_y = y_index - low_y;
	const real frac_z = z_index - low_z;

	return (1 - frac_x)*(1 - frac_y)*(1 - frac_z)*v000
		+ frac_x*(1 - frac_y)*(1 - frac_z)*v100
		+ (1 - frac_x)*frac_y*(1 - frac_z)*v010
		+ frac_x*frac_y*(1 - frac_z)*v110
		+ (1 - frac_x)*(1 - frac_y)*frac_z*v001
		+ frac_x*(1 - frac_y)*frac_z*v101
		+ (1 - frac_x)*frac_y*frac_z*v011
		+ frac_x*frac_y*frac_z*v111;
}

template<typename T, bool gpu_flag>
PHYSICS T table_3D<T, gpu_flag>::get_rounddown(real x, real y, real z) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;
	const real z_index = (z - _z_min) * _z_step;
	const int low_x = clampi(floorr(x_index), 0, _width - 1);
	const int low_y = clampi(floorr(y_index), 0, _height - 1);
	const int low_z = clampi(floorr(z_index), 0, _depth - 1);
	return (*this)(low_x, low_y, low_z);
}

template<typename T, bool gpu_flag>
PHYSICS T table_3D<T, gpu_flag>::get_nearest(real x, real y, real z) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;
	const real z_index = (z - _z_min) * _z_step;
	const int near_x = clampi(rintr(x_index), 0, _width - 1);
	const int near_y = clampi(rintr(y_index), 0, _height - 1);
	const int near_z = clampi(rintr(z_index), 0, _depth - 1);
	return (*this)(near_x, near_y, near_z);
}

template<typename T>
struct _table_3D_factory<T, false>
{
	inline static CPU void allocate(table_3D<T, false> & table, size_t width, size_t height, size_t depth)
	{
		table._width = width;
		table._height = height;
		table._depth = depth;
		table._data = new T[width * height * depth];
		table._pitch = width * sizeof(T);
	}

	inline static CPU void set(table_3D<T, false> & table, T* data)
	{
		memcpy(table._data, data, table._width * table._height * table._depth * sizeof(T));
	}

	template<typename callback_function>
	inline static CPU void mem_scope(table_3D<T, false> & table, callback_function callback)
	{
		// Make indirect arrays
		T** host_pp = new T*[table._height * table._depth];
		for (size_t y = 0; y < table._height*table._depth; ++y)
			host_pp[y] = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(table._data) + y * table._pitch);
		T*** host_ppp = new T**[table._depth];
		for (size_t z = 0; z < table._depth; ++z)
			host_ppp[z] = &host_pp[z * table._height];

		callback(host_ppp);

		delete[] host_ppp;
		delete[] host_pp;
	}

	inline static CPU void free(table_3D<T, false> & table)
	{
		delete[] table._data;
		table._data = nullptr;
		table._width = 0;
		table._height = 0;
		table._depth = 0;
		table._pitch = 0;
	}
};

#if CUDA_AVAILABLE
template<typename T>
struct _table_3D_factory<T, true>
{
	inline static CPU void allocate(table_3D<T, true> & table, size_t width, size_t height, size_t depth)
	{
		table._width = width;
		table._height = height;
		table._depth = depth;
		cuda::cuda_new_3D(&table._data, &table._pitch, width, height, depth);
	}

	inline static CPU void set(table_3D<T, true> & table, T* data)
	{
		const auto width = table._width;
		const auto height = table._height;
		const auto depth = table._depth;
		cuda::cuda_mem_scope_3D(table._data, table._pitch, width, height, depth,
			[data, width, height, depth](T*** device)
			{
				for (size_t z = 0; z < depth; ++z)
					for (size_t y = 0; y < height; ++y)
						for (size_t x = 0; x < width; ++x)
							device[z][y][x] = data[z*width*height + y*width + x];
			});
	}

	template<typename callback_function>
	inline static CPU void mem_scope(table_3D<T, true> & table, callback_function callback)
	{
		cuda::cuda_mem_scope_3D<T>(table._data, table._pitch, table._width, table._height, table._depth, callback);
	}

	inline static CPU void free(table_3D<T, true> & table)
	{
		cudaFree(table._data);
		table._data = nullptr;
		table._width = 0;
		table._height = 0;
		table._depth = 0;
		table._pitch = 0;
	}
};
#endif // CUDA_AVAILABLE

}} // namespace nbl::util

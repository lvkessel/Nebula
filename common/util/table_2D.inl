#include <cstring>

namespace nbl { namespace util {

template<typename T, bool gpu_flag>
HOST table_2D<T, gpu_flag> table_2D<T, gpu_flag>::create(T* data,
	real x_min, real x_max, size_t width,
	real y_min, real y_max, size_t height)
{
	table_2D<T, gpu_flag> table;
	_table_2D_factory<T, gpu_flag>::allocate(table, width, height);
	_table_2D_factory<T, gpu_flag>::set(table, data);

	table._x_min = x_min;
	table._x_step = (width - 1) / (x_max - x_min);
	table._y_min = y_min;
	table._y_step = (height - 1) / (y_max - y_min);

	return table;
}
template<typename T, bool gpu_flag>
HOST void table_2D<T, gpu_flag>::destroy(table_2D<T, gpu_flag> & table)
{
	_table_2D_factory<T, gpu_flag>::free(table);
}

template<typename T, bool gpu_flag>
PHYSICS T & table_2D<T, gpu_flag>::operator()(size_t i, size_t j)
{
	return *(reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(_data) + j * _pitch) + i);
}

template<typename T, bool gpu_flag>
PHYSICS T const & table_2D<T, gpu_flag>::operator()(size_t i, size_t j) const
{
	return *(reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(_data) + j * _pitch) + i);
}

template<typename T, bool gpu_flag>
PHYSICS T table_2D<T, gpu_flag>::get(real x, real y) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;

	const int low_x = clampi(floorr(x_index), 0, _width - 2);
	const int low_y = clampi(floorr(y_index), 0, _height - 2);

	const T v00 = (*this)(low_x, low_y);
	const T v10 = (*this)(low_x + 1, low_y);
	const T v01 = (*this)(low_x, low_y + 1);
	const T v11 = (*this)(low_x + 1, low_y + 1);

	const real frac_x = x_index - low_x;
	const real frac_y = y_index - low_y;

	return (1 - frac_x)*(1 - frac_y)*v00
		+ frac_x*(1 - frac_y)*v10
		+ (1 - frac_x)*frac_y*v01
		+ frac_x*frac_y*v11;
}

template<typename T, bool gpu_flag>
PHYSICS T table_2D<T, gpu_flag>::get_rounddown(real x, real y) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;
	const int low_x = clampi(floorr(x_index), 0, _width - 1);
	const int low_y = clampi(floorr(y_index), 0, _height - 1);
	return (*this)(low_x, low_y);
}

template<typename T, bool gpu_flag>
PHYSICS T table_2D<T, gpu_flag>::get_nearest(real x, real y) const
{
	const real x_index = (x - _x_min) * _x_step;
	const real y_index = (y - _y_min) * _y_step;
	const int near_x = clampi(rintr(x_index), 0, _width - 1);
	const int near_y = clampi(rintr(y_index), 0, _height - 1);
	return (*this)(near_x, near_y);
}

template<typename T>
struct _table_2D_factory<T, false>
{
	inline static HOST void allocate(table_2D<T, false> & table, size_t width, size_t height)
	{
		table._width = width;
		table._height = height;
		table._data = new T[width * height];
		table._pitch = width * sizeof(T);
	}

	inline static HOST void set(table_2D<T, false> & table, T* data)
	{
		memcpy(table._data, data, table._width * table._height * sizeof(T));
	}

	inline static HOST void free(table_2D<T, false> & table)
	{
		delete[] table._data;
		table._data = nullptr;
		table._width = 0;
		table._height = 0;
		table._pitch = 0;
	}
};

#if CUDA_AVAILABLE
template<typename T>
struct _table_2D_factory<T, true>
{
	inline static HOST void allocate(table_2D<T, true> & table, size_t width, size_t height)
	{
		table._width = width;
		table._height = height;
		cuda::cuda_new_2D(&table._data, &table._pitch, width, height);
	}

	inline static HOST void set(table_2D<T, true> & table, T* data)
	{
		const auto width = table._width;
		const auto height = table._height;
		cuda::cuda_mem_scope_2D(table._data, table._pitch, width, height,
			[data, width, height](T** device)
			{
				for (size_t y = 0; y < height; ++y)
					for (size_t x = 0; x < width; ++x)
						device[y][x] = data[y*width + x];
			});
	}

	inline static HOST void free(table_2D<T, true> & table)
	{
		cudaFree(table._data);
		table._data = nullptr;
		table._width = 0;
		table._height = 0;
		table._pitch = 0;
	}
};
#endif // CUDA_AVAILABLE

}} // namespace nbl::util

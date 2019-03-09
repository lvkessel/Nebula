#include <cstring>

namespace nbl { namespace util {

template<typename T, bool gpu_flag>
CPU table_1D<T, gpu_flag> table_1D<T, gpu_flag>::create(
	real x_min, real x_max, size_t n,
	T* data)
{
	table_1D<T, gpu_flag> table;
	detail::table_1D_factory<T, gpu_flag>::allocate(table, n);
	if (data != nullptr)
		detail::table_1D_factory<T, gpu_flag>::set(table, data);
	table._x_min = x_min;
	table._x_step = (n - 1) / (x_max - x_min);
	return table;
}
template<typename T, bool gpu_flag>
CPU void table_1D<T, gpu_flag>::destroy(table_1D<T, gpu_flag> & table)
{
	detail::table_1D_factory<T, gpu_flag>::free(table);
}

template<typename T, bool gpu_flag>
template<typename callback_function>
CPU void table_1D<T, gpu_flag>::mem_scope(callback_function callback)
{
	detail::table_1D_factory<T, gpu_flag>::mem_scope(*this, callback);
}

template<typename T, bool gpu_flag>
PHYSICS T & table_1D<T, gpu_flag>::operator()(size_t i)
{
	return _data[i];
}
template<typename T, bool gpu_flag>
PHYSICS T const & table_1D<T, gpu_flag>::operator()(size_t i) const
{
	return _data[i];
}

template<typename T, bool gpu_flag>
PHYSICS T table_1D<T, gpu_flag>::get(real x) const
{
	const real x_index = (x - _x_min) * _x_step;

	const int low_index = clampi(floorr(x_index), 0, _n - 2);
	const T low_value = _data[low_index];
	const T high_value = _data[low_index + 1];

	const real frac_index = x_index - low_index;
	return (1 - frac_index)*low_value + frac_index*high_value;

/*
Original implementation does not handle infinities correctly.
	return low_value + (x_index - low_index) * (high_value - low_value);
*/
}

namespace detail
{
	template<typename T>
	struct table_1D_factory<T, false>
	{
		inline static CPU void allocate(table_1D<T, false> & table, size_t n)
		{
			table._n = n;
			table._data = new T[n];
		}

		inline static CPU void set(table_1D<T, false> & table, T* data)
		{
			memcpy(table._data, data, table._n * sizeof(T));
		}

		template<typename callback_function>
		inline static CPU void mem_scope(table_1D<T, false> & table, callback_function callback)
		{
			callback(table._data);
		}

		inline static CPU void free(table_1D<T, false> & table)
		{
			delete[] table._data;
			table._data = nullptr;
			table._n = 0;
		}
	};

#if CUDA_COMPILER_AVAILABLE
	template<typename T>
	struct table_1D_factory<T, true>
	{
		inline static CPU void allocate(table_1D<T, true> & table, size_t n)
		{
			table._n = n;
			cuda::cuda_new<T>(&table._data, n);
		}

		inline static CPU void set(table_1D<T, true> & table, T* data)
		{
			const auto n = table._n;
			cuda::cuda_mem_scope<T>(table._data, table._n, [data, n](T* device)
			{
				for (size_t i = 0; i < n; ++i)
					device[i] = data[i];
			});
		}

		template<typename callback_function>
		inline static CPU void mem_scope(table_1D<T, true> & table, callback_function callback)
		{
			cuda::cuda_mem_scope<T>(table._data, table._n, callback);
		}

		inline static CPU void free(table_1D<T, true> & table)
		{
			cudaFree(table._data);
			table._data = nullptr;
			table._n = 0;
		}
	};
#endif // CUDA_COMPILER_AVAILABLE
} // namespace detail

}} // namespace nbl::util

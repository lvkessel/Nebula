#ifndef __TABLE_1D_H_
#define __TABLE_1D_H_

namespace nbl { namespace util {

namespace detail
{
	/**
	 * \brief Responsible for managing the memory on the CPU or GPU device.
	 *
	 * Reason for putting this in a separate struct is that we can't
	 * partially specialize member functions of the table_1D class.
	 */
	template<typename T, bool gpu_flag>
	struct table_1D_factory;
}

/**
 * \brief 1D table to be used in simulation code.
 *
 * It doesn't just store data, but also information about the "axis".
 *
 * The table is allocated by the static {@link create} function and deallocated
 * by the {@link destroy} function. It is up to the user to call these functions
 * correctly! The reason is that these tables are expected to be passed to the
 * GPU as parameters, so we don't want to automatically call constructors and
 * destructors.
 *
 * \tparam T        Datatype stored (usually ::real)
 * \tparam gpu_flag Set to true if data should be stored on the GPU, false for CPU.
 */
template<typename T, bool gpu_flag>
class table_1D
{
public:
	/**
	 * \brief Create a table, allocating memory for the data.
	 *
	 * Throws `std::bad_alloc` exception if memory allocation fails.
	 *
	 * \param x_min Lower value of the physical range
	 * \param x_max Upper value of the physical range
	 * \param n     Number of data points
	 * \param data  Optionally: provide data.
	 */
	static CPU table_1D<T, gpu_flag> create(
		real x_min, real x_max, size_t n,
		T* data = nullptr);

	/**
	 * \brief Destroy a table, deallocating the data.
	 */
	static CPU void destroy(table_1D<T, gpu_flag> & arr);

	/**
	 * \brief Set or get data using a callback function.
	 *
	 * Similar to ::nbl::cuda::cuda_mem_scope. If applicable, it copies the
	 * data from the GPU to the CPU. It calls the \p callback function with a
	 * CPU-accessible pointer to this data. Finally, if applicable, it copies
	 * the data back to the GPU.
	 *
	 * \param callback Callback function. Should have signature `void(T*)`.
	 */
	template<typename callback_function>
	CPU void mem_scope(callback_function callback);

	/**
	 * \brief Direct read-write access to data. No bounds checking is done.
	 *
	 * \param i Index to be accessed.
	 */
	inline PHYSICS T & operator()(size_t i);
	/**
	 * \brief Direct read-only access to data. No bounds checking is done.
	 *
	 * \param i Index to be accessed.
	 */
	inline PHYSICS T const & operator()(size_t i) const;

	/**
	 * \brief Get value at some x coordinate, with linear interpolation and
	 * extrapolation.
	 *
	 * If \p x is between `x_min` and `x_max` provided to the constructor, this
	 * function performs linear interpolation. If it is outside that range, this
	 * function performs linear extrapolation.
	 *
	 * \param x Position to get the data for.
	 */
	inline PHYSICS T get(real x) const;

private:
	T* _data = nullptr;

	size_t _n = 0;
	real _x_min = 0;
	real _x_step = 0;

	friend struct detail::table_1D_factory<T, gpu_flag>;
};

}} // namespace nbl::util

#include "table_1D.inl"

#endif // __TABLE_1D_H_

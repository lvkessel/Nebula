#ifndef __TABLE_1D_H_
#define __TABLE_1D_H_

/*
 * 1D table to be used in simulation code.
 * 
 * The table is allocated by the static create() function and
 * deallocated by the destroy() function. It is up to the user
 * to call these functions correctly!
 * 
 * Template parameters are the datatype stored (usually real)
 * and a flag indicating if the data should be stored on the GPU.
 */

namespace nbl { namespace util {

template<typename T, bool gpu_flag>
struct _table_1D_factory;

template<typename T, bool gpu_flag>
class table_1D
{
public:
	// Throws bad_alloc exception if it goes wrong
	static HOST table_1D<T, gpu_flag> create(T* data, real x_min, real x_max, size_t n);
	static HOST void destroy(table_1D<T, gpu_flag> & arr);

	// Direct access to index, no bounds checking.
	inline PHYSICS T & operator()(size_t i);
	inline PHYSICS T const & operator()(size_t i) const;

	// Get value at some x coordinate, with linear interpolation and extrapolation
	inline PHYSICS T get(real x) const;

private:
	T* _data = nullptr;

	size_t _n = 0;
	real _x_min = 0;
	real _x_step = 0;

	friend struct _table_1D_factory<T, gpu_flag>;
};

}} // namespace nbl::util

#include "table_1D.inl"

#endif // __TABLE_1D_H_

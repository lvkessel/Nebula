#ifndef __TABLE_2D_H_
#define __TABLE_2D_H_

/*
 * 2D table to be used in simulation code.
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
struct _table_2D_factory;

template<typename T, bool gpu_flag>
class table_2D
{
public:
	// Data organised such that data[x,y] = data + y*width + x
	// Throws bad_alloc exception if it goes wrong
	static HOST table_2D<T, gpu_flag> create(T* data,
		real x_min, real x_max, size_t width,
		real y_min, real y_max, size_t height);
	static HOST void destroy(table_2D<T, gpu_flag> & arr);

	// Direct access to index, no bounds checking.
	inline PHYSICS T & operator()(size_t i, size_t j);
	inline PHYSICS T const & operator()(size_t i, size_t j) const;

	// Get value at some x coordinate, with linear interpolation and extrapolation
	inline PHYSICS T get(real x, real y) const;
	// Convert position to indices, round them down, and return non-interpolated value.
	inline PHYSICS T get_rounddown(real x, real y) const;
	// Convert position to indices, round to nearest, and return non-interpolated value.
	inline PHYSICS T get_nearest(real x, real y) const;

private:
	T* _data = nullptr;
	size_t _pitch = 0;

	size_t _width = 0;
	real _x_min = 0;
	real _x_step = 0;

	size_t _height = 0;
	real _y_min = 0;
	real _y_step = 0;

	friend struct _table_2D_factory<T, gpu_flag>;
};

}} // namespace nbl::util

#include "table_2D.inl"

#endif // __TABLE_2D_H_

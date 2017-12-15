#ifndef __RANDOM_H_
#define __RANDOM_H_

/*
 * The _random_state class holds the random number generator, which is CPU or
 * GPU specific.
 * Users should create a random_generator class, with template parameter
 * gpu_flag ==  true for a GPU random number generator or false for a CPU random
 * number generator. It inherits from _random_state to provide device-specific
 * functionality.
 */

namespace nbl { namespace util {

template<bool gpu_flag>
class _random_state;

template<bool gpu_flag>
class random_generator
	: public _random_state<gpu_flag>
{
public:
	using _random_state<gpu_flag>::_random_state;

	// Uniformly distributed between 0 and 1
	using _random_state<gpu_flag>::unit;

	// Uniformly distributed between 0 and 2*pi
	PHYSICS real phi();

	// Exponential, with typical constant tau (tau has units of return value)
	PHYSICS real exponential(real tau);
};

}} // namespace nbl::random

#include "random.inl"

#endif // __RANDOM_H_

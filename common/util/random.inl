#include <random>

namespace nbl { namespace util {

template<>
class _random_state<false>
{
public:
	HOST _random_state() {}
	PHYSICS real unit()
	{
#if CUDA_COMPILING
		// TODO: proper error message.
		return 0;
#else
		return std::generate_canonical<real, std::numeric_limits<real>::digits>(_generator);
#endif
	}

private:
	std::mt19937 _generator;
};

template<bool gpu_flag>
PHYSICS real random_generator<gpu_flag>::phi()
{
	return 2 * pi * unit();
}

template<bool gpu_flag>
PHYSICS real random_generator<gpu_flag>::exponential(real tau)
{
	return -tau*logr(unit());
}


#if CUDA_AVAILABLE
#include "curand_kernel.h"

template<>
class _random_state<true>
{
public:
	__device__ _random_state(unsigned long long seed, unsigned long long sequence)
	{
		curand_init(seed, sequence, 0, &_rand_state);
	}
	PHYSICS real unit()
	{
#if CUDA_COMPILING
#if USE_DOUBLE
		return curand_uniform_double(&_rand_state);
#else // USE_DOUBLE
		return curand_uniform(&_rand_state);
#endif // USE_DOUBLE
#else // CUDA_COMPILING
		// TODO: proper error message.
		return 0;
#endif // CUDA_COMPILING
	}

private:
	curandState _rand_state;
};

#endif // CUDA_AVAILABLE

}} // namespace nbl::util
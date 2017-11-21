#ifndef __CUDA_NEW_H_
#define __CUDA_NEW_H_

namespace nbl { namespace cuda {

template<typename T>
HOST void cuda_new(T** destination, size_t N)
{
	cudaError_t cudaStatus = cudaMalloc(destination, N * sizeof(T));
	if (cudaStatus != cudaSuccess)
	{
		throw std::bad_alloc();
	}
}

template<typename T>
HOST void cuda_new_2D(T** destination, size_t* pitch, size_t width, size_t height)
{
	cudaError_t cudaStatus = cudaMallocPitch(destination, pitch, width*sizeof(T), height);
	if (cudaStatus != cudaSuccess)
	{
		throw std::bad_alloc();
	}
}

template<typename T>
HOST void cuda_new_3D(T** destination, size_t* pitch, size_t width, size_t height, size_t depth)
{
	cudaError_t cudaStatus = cudaMallocPitch(destination, pitch, width*sizeof(T), height*depth);
	if (cudaStatus != cudaSuccess)
	{
		throw std::bad_alloc();
	}
}

}} // namespace nbl::cuda

#endif // __CUDA_NEW_H_
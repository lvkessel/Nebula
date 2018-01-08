#ifndef __CUDA_MEM_SCOPE_H_
#define __CUDA_MEM_SCOPE_H_

namespace nbl { namespace cuda {

/*
 * Read/write data in GPU memory: copy to host, callback sets data, copy to device.
 * Callback function should have signature void(T*).
 */
template<typename T, typename callback_func>
CPU void cuda_mem_scope(T* dev_p, size_t N, callback_func callback)
{
	if (dev_p == nullptr || N == 0)
		return;

	// Copy to host.
	// Use malloc, not new, to avoid calling constructors.
	T* host_p = reinterpret_cast<T*>(malloc(N * sizeof(T)));
	cudaMemcpy(host_p, dev_p, N * sizeof(T), cudaMemcpyDeviceToHost);

	// Callback
	callback(host_p);

	// Copy back to device
	cudaMemcpy(dev_p, host_p, N * sizeof(T), cudaMemcpyHostToDevice);
	free(host_p);
}

/*
 * Similar to cuda_mem_scope, only 2D access this time.
 * Callback function should have signature void(T** arr),
 * with indexing convention arr[y][x]
 */
template<typename T, typename callback_func>
CPU void cuda_mem_scope_2D(T* dev_p, size_t pitch, size_t width, size_t height, callback_func callback)
{
	if (dev_p == nullptr || pitch == 0 || width == 0 || height == 0)
		return;

	// Copy to host
	T* host_p = reinterpret_cast<T*>(new uint8_t[pitch * height]);
	cudaMemcpy(host_p, dev_p, pitch*height, cudaMemcpyDeviceToHost);

	// Make indirect array
	T** host_pp = new T*[height];
	for (size_t y = 0; y < height; ++y)
		host_pp[y] = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(host_p) + y*pitch);

	// Callback
	callback(host_pp);

	// Copy back
	cudaMemcpy(dev_p, host_p, pitch*height, cudaMemcpyHostToDevice);
	
	delete[] host_pp;
	delete[] host_p;
}

/*
 * Same, for 3D.
 * Callback function should have signature void(T*** arr),
 * with indexing convention arr[z][y][x]
 */
template<typename T, typename callback_func>
CPU void cuda_mem_scope_3D(T* dev_p, size_t pitch, size_t width, size_t height, size_t depth, callback_func callback)
{
	if (dev_p == nullptr || pitch == 0 || width == 0 || height == 0 || depth == 0)
		return;

	// Copy to host
	T* host_p = reinterpret_cast<T*>(new uint8_t[pitch * height * depth]);
	cudaMemcpy(host_p, dev_p, pitch*height*depth, cudaMemcpyDeviceToHost);

	// Make indirect arrays
	T** host_pp = new T*[height*depth];
	for (size_t y = 0; y < height*depth; ++y)
		host_pp[y] = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(host_p) + y*pitch);
	T*** host_ppp = new T**[depth];
	for (size_t z = 0; z < depth; ++z)
		host_ppp[z] = host_pp[z * height];

	// Callback
	callback(host_ppp);

	// Copy back
	cudaMemcpy(dev_p, host_p, pitch*height*depth, cudaMemcpyHostToDevice);

	delete[] host_ppp;
	delete[] host_pp;
	delete[] host_p;
}

}} // namespace nbl::cuda

#endif // __CUDA_MEM_SCOPE_H_

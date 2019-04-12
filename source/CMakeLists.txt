add_subdirectory(legacy_thomas)
add_subdirectory(geometry/octree)

if(CMAKE_CUDA_COMPILER)
	add_executable(nebula_gpu gpu_main.cu material/hdf5_file.cpp
		$<TARGET_OBJECTS:legacy_thomas>
		$<TARGET_OBJECTS:octree>)
	target_link_libraries(nebula_gpu ${LIBRARIES})
endif()

add_executable(nebula_cpu_mt cpu_mt_main.cpp material/hdf5_file.cpp
		$<TARGET_OBJECTS:legacy_thomas>
		$<TARGET_OBJECTS:octree>)
add_executable(nebula_cpu_inspect cpu_inspect_main.cpp material/hdf5_file.cpp
		$<TARGET_OBJECTS:legacy_thomas>
		$<TARGET_OBJECTS:octree>)
target_link_libraries(nebula_cpu_mt ${LIBRARIES} Threads::Threads)
target_link_libraries(nebula_cpu_inspect ${LIBRARIES})

add_subdirectory(legacy_thomas)

add_executable(nebula_gpu gpu_main.cu material/hdf5_file.cpp)
add_executable(nebula_cpu_mt cpu_mt_main.cpp material/hdf5_file.cpp)
add_executable(nebula_cpu_inspect cpu_inspect_main.cpp material/hdf5_file.cpp)

target_link_libraries(nebula_gpu legacy_thomas ${LIBRARIES})
target_link_libraries(nebula_cpu_mt legacy_thomas ${LIBRARIES} Threads::Threads)
target_link_libraries(nebula_cpu_inspect legacy_thomas ${LIBRARIES})

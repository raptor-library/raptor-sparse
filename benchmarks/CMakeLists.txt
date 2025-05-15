include_directories(${raptor_sparse_INCDIR})

add_executable(benchmark_spmv benchmark_spmv.cpp)
target_link_libraries(benchmark_spmv raptor-sparse)

if (WITH_MPI)
    add_executable(benchmark_par_spmv benchmark_par_spmv.cpp)
    target_link_libraries(benchmark_par_spmv raptor-sparse ${MPI_LIBRARIES})

    add_executable(benchmark_par_block_spmv benchmark_par_block_spmv.cpp)
    target_link_libraries(benchmark_par_block_spmv raptor-sparse ${MPI_LIBRARIES})
endif()



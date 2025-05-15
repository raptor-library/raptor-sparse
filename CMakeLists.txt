cmake_minimum_required(VERSION 3.12...3.20)
project(raptor-sparse)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wredundant-decls -Wcast-align -Wshadow")

add_compile_definitions(RAPTOR_SPARSE_TEST_FOLDER="${CMAKE_SOURCE_DIR}/src/tests/")

include(FeatureSummary)

option(WITH_PTSCOTCH "Enable PTScotch Partitioning" OFF)
option(WITH_PARMETIS "Enable ParMetis Partitioning" OFF)
option(WITH_MFEM "Add MFEM" OFF)
option(WITH_MPI "Using MPI" ON)
option(WITH_HOSTFILE "Use a Hostfile with MPI" OFF)

add_feature_info(mfem WITH_MFEM "MFEM matrix gallery")
add_feature_info(ptscotch WITH_PTSCOTCH "Enable PTScotch Partitioning")
add_feature_info(parmetis WITH_PARMETIS "Enable ParMetis Partitioning")
add_feature_info(hostfile WITH_HOSTFILE "Enable Hostfile for MPIRUN")

include(options)
include(testing)
include(source)
SET(MPIRUN "")

if (WITH_MPI)
    add_definitions ( -DUSING_MPI )
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    SET(MPIRUN mpirun)
endif (WITH_MPI)

set(raptor_sparse_DIR src)
set(raptor_sparse_INCDIR ${CMAKE_CURRENT_SOURCE_DIR}/${raptor_sparse_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

SET(MFEM_DIR "" CACHE STRING "Directory of mfem")
SET(MFEM_MESH_DIR "" CACHE STRING "Directory containing MFEM mesh files")
SET(PTSCOTCH_DIR "" CACHE STRING "Directory of Ptscotch")
SET(PARMETIS_DIR "" CACHE STRING "Directory of ParMetis")
SET(HOSTFILE "" CACHE STRING "Set Hostfile")
SET(BLAS_PATH "" CACHE STRING "/usr/linb/x86_64-linux-gnu")
SET(LAPACK_PATH "" CACHE STRING "/usr/linb/x86_64-linux-gnu")

find_library(LAPACK_LIB NAMES liblapack.so.3 lapack HINTS ${LAPACK_PATH})
find_library(BLAS_LIB NAMES libblas.so.3 blas HINTS ${BLAS_PATH})
set(EXTERNAL_LIBS ${LAPACK_LIB} ${BLAS_LIB})

if (WITH_HOSTFILE)
    find_file (FILE_OF_HOST, ${HOSTFILE})
    set(HOST "--hostfile" "${HOSTFILE}")
endif()

if (WITH_PTSCOTCH)
    add_definitions ( -DUSING_PTSCOTCH )
    find_package(Ptscotch)
    if (PTSCOTCH_FOUND)
        include_directories(${PTSCOTCH_INCLUDE_DIRS})
    else()
        message(FATAL_ERROR "Cannot find PTSCOTCH.  Try setting PTSCOTCH_DIR" )
        set(WITH_PTSCOTCH ON)
    endif(PTSCOTCH_FOUND)
endif(WITH_PTSCOTCH)

if (WITH_PARMETIS)
    add_definitions ( -DUSING_PARMETIS )
    find_package(ParMetis)
    if (PARMETIS_FOUND)
        include_directories(${PARMETIS_INCLUDE_DIRS})
    else()
        message(FATAL_ERROR "Cannot find PARMETIS.  Try setting PARMETIS_DIR" )
        set(WITH_PARMETIS ON)
    endif(PARMETIS_FOUND)
endif(WITH_PARMETIS)

if(WITH_MFEM)
    add_definitions ( -DUSING_MFEM )
	find_package(Mfem)
	if (MFEM_FOUND)
		include_directories(${MFEM_INCLUDE_DIRS})
		message(STATUS "MFEM Mesh Directory: " ${MFEM_MESH_PATH})
	else()
		message(FATAL_ERROR "Cannot find mfem. Try setting MFEM_DIR.")
	endif(MFEM_FOUND)
    add_definitions( -DMFEM_MESH_DIR=\"${MFEM_MESH_PATH}\" )
endif(WITH_MFEM)

add_subdirectory(src)

if (BUILD_EXAMPLES)
  add_subdirectory(benchmarks)
endif()

find_package(Doxygen)

if (DOXYGEN_FOUND)
    # make the output directory
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)

    # copy the doxygen file with @whatever@ substitutions
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "generating Doxygen output"
        VERBATIM
        )
endif (DOXYGEN_FOUND)

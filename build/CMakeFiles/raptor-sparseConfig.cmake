

if(NOT TARGET raptor-sparse AND NOT raptor-sparse_BINARY_DIR)
  include("${CMAKE_CURRENT_LIST_DIR}/raptor-sparseTargets.cmake")
  endif()

find_package(Threads REQUIRED)
find_package(MPI COMPONENTS CXX REQUIRED)

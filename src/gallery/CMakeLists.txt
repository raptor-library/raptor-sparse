# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if (WITH_MFEM)
    set(mfem_HEADERS
        external/mfem_wrapper.hpp
    )
    set(mfem_SOURCES 
        external/mfem/mfem_linear_elasticity.cpp
        external/mfem/mfem_laplacian.cpp
        external/mfem/mfem_grad_div.cpp
        external/mfem/mfem_adaptive_laplacian.cpp
        external/mfem/mfem_dg_diffusion.cpp
        external/mfem/mfem_dg_elasticity.cpp
       )
endif()

# Create a variable called gallery_SOURCES containing all .cpp files:
if (WITH_MPI)
    set(par_gallery_HEADERS
        gallery/par_stencil.hpp
        gallery/par_random.hpp
        gallery/par_matrix_IO.hpp
        gallery/par_matrix_market.hpp
        )
    set(par_gallery_SOURCES
        gallery/par_stencil.cpp
        gallery/par_random.cpp
        gallery/par_matrix_IO.cpp
        gallery/par_matrix_market.cpp
        )
else ()
    set(par_gallery_HEADERS
        ""
        )
    set(par_gallery_SOURCES
        ""
        )
endif()

set(gallery_HEADERS
    gallery/diffusion.hpp
    gallery/laplacian27pt.hpp
    gallery/stencil.hpp
    gallery/random.hpp
    gallery/matrix_IO.hpp
    gallery/matrix_market.hpp
    ${par_gallery_HEADERS}
    PARENT_SCOPE)

set(gallery_SOURCES
    gallery/diffusion.cpp
    gallery/laplacian27pt.cpp
    gallery/stencil.cpp
    gallery/random.cpp
    gallery/matrix_IO.cpp
    gallery/matrix_market.cpp
    ${par_gallery_SOURCES}
    PARENT_SCOPE)



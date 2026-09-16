set(curdir ${CMAKE_CURRENT_SOURCE_DIR})

set(sources-fpp
  ${curdir}/arpack.F90
  ${curdir}/blas.F90
  ${curdir}/chimes.F90
  ${curdir}/ddcosmo.F90
  ${curdir}/dftd4refs.F90
  ${curdir}/elsiiface.F90
  ${curdir}/lapack.F90
  ${curdir}/lebedev.F90
  ${curdir}/magma.F90
  ${curdir}/mpifx.F90
  ${curdir}/plumed.F90
  ${curdir}/poisson.F90
  ${curdir}/scalapackfx.F90
  ${curdir}/sdftd3.F90
  ${curdir}/tblite.F90
  ${curdir}/xmlf90.F90)

if(WITH_SOCKETS)
  list(APPEND sources-fpp ${curdir}/fsockets.F90)
endif(WITH_SOCKETS)

if(WITH_TRANSPORT)
  list(APPEND sources-fpp
    ${curdir}/negf.F90)
endif()

if(WITH_MAGMA)
  list(APPEND sources-c
    ${curdir}/magmac.c)
endif()

if(WITH_MBD)
  list(APPEND sources-fpp ${curdir}/mbd.F90)
endif()

# Ewald engines copied from ewald_summation; shared modules are compiled once.
set(sources-fftewald-f90
  ${curdir}/fftewald/common/ewald_alpha_search.f90
  ${curdir}/fftewald/common/ewald_cell_list.f90
  ${curdir}/fftewald/common/ewald_constants.f90
  ${curdir}/fftewald/common/ewald_geometry.f90
  ${curdir}/fftewald/common/ewald_lambert_w.f90
  ${curdir}/fftewald/common/ewald_self.f90
  ${curdir}/fftewald/common/ewald_truncation.f90
  ${curdir}/fftewald/common/ewald_validation.f90
  ${curdir}/fftewald/direct/ewald_direct_2d.f90
  ${curdir}/fftewald/direct/ewald_direct_3d.f90
  ${curdir}/fftewald/fft/ewald_fft_2d.f90
  ${curdir}/fftewald/fft/ewald_fft_2d_fourier.f90
  ${curdir}/fftewald/fft/ewald_fft_2d_parameters.f90
  ${curdir}/fftewald/fft/ewald_fft_2d_real.f90
  ${curdir}/fftewald/fft/ewald_fft_3d.f90
  ${curdir}/fftewald/fft/ewald_fft_3d_fourier.f90
  ${curdir}/fftewald/fft/ewald_fft_3d_parameters.f90
  ${curdir}/fftewald/fft/ewald_fft_3d_real.f90
  ${curdir}/fftewald/fft/ewald_slab_kernel.f90
  ${curdir}/fftewald/nfft/fft_backend.f90
  ${curdir}/fftewald/nfft/nfft.f90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
set(ALL-SOURCES-F90 ${ALL-SOURCES-F90} ${sources-fftewald-f90} PARENT_SCOPE)
set(ALL-SOURCES-C ${ALL-SOURCES-C} ${sources-c} PARENT_SCOPE)

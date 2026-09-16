module fft_backend
   !> Uniform fast Fourier transforms, provided by FFTW3.  This is the
   !> library's only external dependency, and isolating it here means the rest
   !> of the code never mentions FFTW or handles a plan.
   !>
   !> Sign convention: the transforms take an integer whose sign selects the
   !> direction, non-negative for the sign +1 convention used by the adjoint
   !> non-uniform transform, negative for the sign -1 of an ordinary forward
   !> transform.  No normalisation is applied in either direction.
   !>
   !> Storage order: FFTW stores arrays with the last index contiguous, Fortran
   !> with the first, so the multi-dimensional planners are given their
   !> dimensions in reverse.  A multi-dimensional transform is separable, so the
   !> per-axis frequency layout is that of transforming one axis at a time,
   !> which is what the deconvolution tables in the nfft module assume.
   use ewald_constants, only: dp
   use, intrinsic :: iso_c_binding
   !$ use omp_lib, only: omp_get_max_threads
   implicit none
   include 'fftw3.f03'

   private
   public :: fft_1d_inplace, fft_2d_inplace, fft_3d_inplace
   public :: fft_2d_real_to_complex, fft_3d_real_to_complex
   public :: next_fast_length, minTransformLength

   !> Shortest transform the library will use along any axis.  Below this the
   !> mode set is too coarse to represent anything and the transform is not
   !> worth its bookkeeping.
   integer, parameter :: minTransformLength = 4

   !> Whether FFTW's threading layer has been initialised.  Initialising it
   !> more than once is an error, so the first transform that needs threads
   !> sets this and every later one sees it.
   logical, save :: threadsInitialised = .false.

contains

   ! =====================================================================
   !  Complex transforms
   ! =====================================================================

   !> One-dimensional in-place transform of a complex vector.
   !>
   !> Deliberately single-threaded and plan-cached.  Its only caller is the
   !> transform of the regularised slab kernel, which is short and is repeated
   !> once per in-plane mode with identical length and direction, so planning
   !> would otherwise dominate the transform itself.  The unaligned flag keeps
   !> the cached plan valid for whatever array the next caller passes.
   !>
   !> The cache makes this routine unsafe to call from several threads at once;
   !> its callers are serial loops.
   subroutine fft_1d_inplace(vector, n, isign)
      integer, intent(in) :: n
      integer, intent(in) :: isign
      complex(dp), intent(inout) :: vector(0:n - 1)

      type(C_PTR), save :: plan = C_NULL_PTR     ! the cached plan
      integer, save :: cachedLength = -1         ! length it was built for
      integer, save :: cachedDirection = 0       ! direction it was built for

      if (n /= cachedLength .or. fftw_direction(isign) /= cachedDirection) then
         if (c_associated(plan)) call fftw_destroy_plan(plan)
         ! The estimating planner is used because it does not overwrite the
         ! input array while planning, which a measuring planner would.
         plan = fftw_plan_dft_1d(int(n, C_INT), vector, vector, &
                                 fftw_direction(isign), &
                                 ior(FFTW_ESTIMATE, FFTW_UNALIGNED))
         cachedLength = n
         cachedDirection = fftw_direction(isign)
      end if

      call fftw_execute_dft(plan, vector, vector)

   end subroutine fft_1d_inplace

   !> Two-dimensional in-place transform of a complex grid.
   subroutine fft_2d_inplace(grid, n1, n2, isign)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: isign
      complex(dp), intent(inout) :: grid(0:n1 - 1, 0:n2 - 1)

      type(C_PTR) :: plan

      call plan_over_all_threads()
      ! Dimensions reversed: see the note on storage order in the module header.
      plan = fftw_plan_dft_2d(int(n2, C_INT), int(n1, C_INT), &
                              grid, grid, fftw_direction(isign), FFTW_ESTIMATE)
      call fftw_execute_dft(plan, grid, grid)
      call fftw_destroy_plan(plan)

   end subroutine fft_2d_inplace

   !> Three-dimensional in-place transform of a complex grid.
   subroutine fft_3d_inplace(grid, n1, n2, n3, isign)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      integer, intent(in) :: isign
      complex(dp), intent(inout) :: grid(0:n1 - 1, 0:n2 - 1, 0:n3 - 1)

      type(C_PTR) :: plan

      call plan_over_all_threads()
      plan = fftw_plan_dft_3d(int(n3, C_INT), int(n2, C_INT), int(n1, C_INT), &
                              grid, grid, fftw_direction(isign), FFTW_ESTIMATE)
      call fftw_execute_dft(plan, grid, grid)
      call fftw_destroy_plan(plan)

   end subroutine fft_3d_inplace

   ! =====================================================================
   !  Real-input transforms
   ! =====================================================================

   !> Two-dimensional in-place transform of a real grid.
   !>
   !> The charge spreading of the non-uniform transform produces a purely real
   !> grid, whose spectrum is Hermitian, so only the indices 0..n1/2 along the
   !> first axis are stored.  That halves the memory the long-range branch
   !> needs.
   !>
   !> The two arrays must alias the same storage for the transform to be in
   !> place: the caller allocates a real grid whose first dimension is padded
   !> to 2*(n1/2+1) and points a complex view at it.
   !>
   !> FFTW's real-input transform is fixed to the sign -1 direction whereas the
   !> adjoint non-uniform transform uses sign +1.  For a real input the two
   !> differ by complex conjugation alone, so the magnitudes an energy needs
   !> are identical; a caller needing the signed structure factor must use the
   !> complex transform instead.
   subroutine fft_2d_real_to_complex(realGrid, complexGrid, n1, n2)
      integer, intent(in) :: n1
      integer, intent(in) :: n2

      !> Padded real grid: logical extent (n1, n2), leading dimension
      !> 2*(n1/2+1).  Overwritten by the transform.
      real(dp), intent(inout) :: realGrid(0:2*(n1/2 + 1) - 1, 0:n2 - 1)

      !> Complex view of the same storage, receiving the half spectrum.
      complex(dp), intent(inout) :: complexGrid(0:n1/2, 0:n2 - 1)

      type(C_PTR) :: plan

      call plan_over_all_threads()
      plan = fftw_plan_dft_r2c_2d(int(n2, C_INT), int(n1, C_INT), &
                                  realGrid, complexGrid, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(plan, realGrid, complexGrid)
      call fftw_destroy_plan(plan)

   end subroutine fft_2d_real_to_complex

   !> Three-dimensional in-place transform of a real grid.  Same idea, same
   !> conventions and the same aliasing requirement as the two-dimensional
   !> version above.
   subroutine fft_3d_real_to_complex(realGrid, complexGrid, n1, n2, n3)
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3

      !> Padded real grid: logical extent (n1, n2, n3), leading dimension
      !> 2*(n1/2+1).  Overwritten by the transform.
      real(dp), intent(inout) :: realGrid(0:2*(n1/2 + 1) - 1, 0:n2 - 1, 0:n3 - 1)

      !> Complex view of the same storage, receiving the half spectrum.
      complex(dp), intent(inout) :: complexGrid(0:n1/2, 0:n2 - 1, 0:n3 - 1)

      type(C_PTR) :: plan

      call plan_over_all_threads()
      plan = fftw_plan_dft_r2c_3d(int(n3, C_INT), int(n2, C_INT), int(n1, C_INT), &
                                  realGrid, complexGrid, FFTW_ESTIMATE)
      call fftw_execute_dft_r2c(plan, realGrid, complexGrid)
      call fftw_destroy_plan(plan)

   end subroutine fft_3d_real_to_complex

   ! =====================================================================
   !  Transform lengths and threading
   ! =====================================================================

   !> Smallest even length not below n whose only prime factors are 2, 3, 5
   !> and 7.
   !>
   !> An FFT is fastest at such lengths, and this library sizes its grids from
   !> a physical cutoff, which lands on no particular integer.  Rounding up
   !> here costs a few percent of grid where rounding up to a power of two
   !> would cost up to a factor of two per axis.  Evenness is required
   !> separately: the mode set runs from -M/2 to M/2-1, which is only well
   !> defined for even M.
   pure function next_fast_length(n) result(m)
      integer, intent(in) :: n
      integer :: m

      integer :: remainder    ! what is left of m after dividing out 2, 3, 5, 7

      m = max(n, 2)
      if (mod(m, 2) /= 0) m = m + 1
      do
         remainder = m
         do while (mod(remainder, 2) == 0); remainder = remainder/2; end do
         do while (mod(remainder, 3) == 0); remainder = remainder/3; end do
         do while (mod(remainder, 5) == 0); remainder = remainder/5; end do
         do while (mod(remainder, 7) == 0); remainder = remainder/7; end do
         if (remainder == 1) exit
         m = m + 2
      end do

   end function next_fast_length

   !> Translate this library's direction convention into FFTW's flag.
   pure function fftw_direction(isign) result(direction)
      integer, intent(in) :: isign
      integer(C_INT) :: direction

      if (isign >= 0) then
         direction = FFTW_BACKWARD
      else
         direction = FFTW_FORWARD
      end if

   end function fftw_direction

   !> Ask FFTW to plan the next transform over as many threads as OpenMP has.
   !> The large grid transforms are always issued from outside a parallel
   !> region, so FFTW's threads never nest inside the library's own.
   subroutine plan_over_all_threads()

      integer :: nThreads

      nThreads = 1
      !$ nThreads = omp_get_max_threads()
      call initialise_threads()
      call fftw_plan_with_nthreads(int(nThreads, C_INT))

   end subroutine plan_over_all_threads

   !> Initialise FFTW's threading layer, exactly once per run.
   subroutine initialise_threads()

      if (.not. threadsInitialised) then
         if (fftw_init_threads() == 0) &
            error stop "fft_backend: FFTW could not initialise its threading layer"
         threadsInitialised = .true.
      end if

   end subroutine initialise_threads

end module fft_backend

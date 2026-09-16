module nfft
   !> Non-uniform fast Fourier transform with a Gaussian window.
   !>
   !> Given N charges q_j at arbitrary positions t_j on the unit torus, the
   !> adjoint transform produces the structure factors
   !>
   !>   Shat(k) = sum_j q_j exp(+2 pi i k . t_j)
   !>
   !> for every integer mode k in a rectangular set; the forward transform is
   !> its exact transpose and evaluates a trigonometric polynomial at those
   !> same positions.  Both cost O(N + M log M) instead of the O(N M) of the
   !> definition, in three steps:
   !>
   !>   1. Spread.  Each charge is smeared onto an oversampled regular grid
   !>      with a Gaussian truncated to a stencil of 2m+1 points per axis.
   !>   2. Transform.  One ordinary FFT of the whole grid.
   !>   3. Deconvolve.  Convolution in real space is multiplication in Fourier
   !>      space, so dividing each retained mode by the Fourier transform of
   !>      the window undoes step 1.  The Gaussian has no zeros, which makes
   !>      that division safe.
   !>
   !> The forward transform runs the same steps backwards.  Accuracy is
   !> controlled by the oversampling factor and the stencil half-width m: the
   !> error falls exponentially in m while the cost grows like m to the power
   !> of the dimension, which is why m is derived from the requested accuracy.
   !>
   !> spread_charges_real_3d is an energy-only shortcut.  The Ewald
   !> energy needs only the squared magnitudes of the structure factors and the
   !> spread grid is purely real, so a caller in that situation can spread onto
   !> a real grid, run a real-input FFT itself and fold the deconvolution into
   !> its own mode sum, which halves the grid memory.  The routine therefore
   !> stops after step 1 and exports the window parameters so that the caller can
   !> finish with exactly the same window.
   !>
   !> Modes come back in the usual FFT bin order: array index a in 0..M-1 holds
   !> the signed mode mode_of_bin(a, M).
   use ewald_constants, only: dp, pi
   use fft_backend, only: fft_3d_inplace
   implicit none

   private
   public :: adjoint_nfft_3d
   public :: forward_nfft_3d
   public :: spread_charges_real_3d
   public :: mode_of_bin, window_shape, window_cutoff_from_budget
   public :: oversampling

   !> Ratio between the working grid and the mode set: the grid has
   !> oversampling*M points per axis for M retained modes.  Two is the standard
   !> choice; every further increase buys less than spending the same memory on
   !> a wider stencil would.
   integer, parameter :: oversampling = 2

   !> Rate at which the window error falls per unit of stencil half-width.  The
   !> truncation error is bounded by a constant times
   !> exp(-pi*m*(1 - 1/(2*oversampling - 1))), which at an oversampling of two
   !> gives 2*pi/3.  Tied to the constant above; recompute it if that changes.
   real(dp), parameter :: windowErrorRate = 2.0_dp*pi/3.0_dp

contains

   ! =====================================================================
   !  Window and index helpers
   ! =====================================================================

   !> Signed mode held by FFT bin a of a transform of length M.
   !>
   !> The bins run 0..M-1 but the modes run -M/2..M/2-1, and an FFT stores the
   !> negative modes in the upper half of the array.
   pure function mode_of_bin(a, M) result(k)

      !> Bin index, in 0..M-1.
      integer, intent(in) :: a

      !> Length of the transform along this axis.
      integer, intent(in) :: M

      !> The signed mode that bin holds.
      integer :: k

      k = a
      if (a >= M/2) k = a - M

   end function mode_of_bin

   !> Shape parameter of the Gaussian window for a stencil of half-width m.
   !> It is chosen so that the window is as narrow as it can be in real space
   !> without its transform decaying too fast across the retained modes: too
   !> wide a window truncates badly, too narrow a one makes the deconvolution
   !> amplify the aliasing error.
   pure function window_shape(m) result(b)

      !> Stencil half-width, in grid points.
      integer, intent(in) :: m

      !> Shape parameter of the Gaussian.
      real(dp) :: b

      b = (2.0_dp*oversampling/(2.0_dp*oversampling - 1.0_dp))*(real(m, dp)/pi)

   end function window_shape

   !> Stencil half-width whose window error sits at or below a given budget.
   !> Solving 4*exp(-windowErrorRate*m) <= budget gives
   !> m = ceil(ln(4/budget)/windowErrorRate), which grows only logarithmically
   !> in the accuracy: six points at 1e-4, fourteen at 1e-12.  The bound is
   !> pessimistic, so the window is never what limits an evaluation.
   pure function window_cutoff_from_budget(budget) result(m)

      !> Share of the requested accuracy the window may spend.
      real(dp), intent(in) :: budget

      !> Stencil half-width, at least two.
      integer :: m

      m = max(2, ceiling(log(4.0_dp/budget)/windowErrorRate))

   end function window_cutoff_from_budget

   !> Gaussian weights of the stencil around one coordinate of one node, and
   !> the index of its lowest grid point.  With the node at t and the grid at
   !> spacing 1/n, the point of index c+o lies at distance (frac - o)/n from
   !> the node, frac being the fractional part of t*n.
   !>
   !> The fractional part must be taken before the grid index is wrapped.  A
   !> fractional coordinate that comes out of a cancelling dot product as a very
   !> small negative number, which happens routinely for a non-orthogonal cell,
   !> wraps to a value just below one; wrapping the index down first would leave
   !> a remainder near n, collapsing every stencil weight to zero and dropping
   !> the charge from the grid.
   subroutine stencil_weights(t, n, m, b, lowestPoint, weights)

      !> Coordinate of the node, interpreted modulo one.
      real(dp), intent(in) :: t

      !> Number of grid points along this axis.
      integer, intent(in) :: n

      !> Stencil half-width in grid points.
      integer, intent(in) :: m

      !> Shape parameter of the Gaussian, from window_shape.
      real(dp), intent(in) :: b

      !> Index of the grid point at offset zero of the stencil.
      integer, intent(out) :: lowestPoint

      !> The 2m+1 window weights, indexed by offset.
      real(dp), intent(out) :: weights(-m:m)

      real(dp) :: wrapped        ! t reduced to [0, 1)
      real(dp) :: scaled         ! wrapped * n, i.e. the node in grid units
      real(dp) :: fraction       ! distance from the node to lowestPoint
      real(dp) :: prefactor      ! normalisation of the Gaussian
      integer  :: offset

      wrapped = t - floor(t)
      scaled = wrapped*real(n, dp)
      lowestPoint = floor(scaled)
      fraction = scaled - real(lowestPoint, dp)     ! in [0, 1) by construction
      if (lowestPoint >= n) lowestPoint = lowestPoint - n   ! wrap only now

      prefactor = 1.0_dp/sqrt(pi*b)
      do offset = -m, m
         weights(offset) = prefactor*exp(-(fraction - real(offset, dp))**2/b)
      end do

   end subroutine stencil_weights

   !> Per-axis tables for the deconvolution step: for each output bin, the fine
   !> grid index it reads from and the factor that undoes the window there.
   subroutine deconvolution_table(M, n, b, gridIndex, factor)

      !> Number of retained modes along this axis.
      integer, intent(in) :: M

      !> Number of fine grid points along this axis.
      integer, intent(in) :: n

      !> Shape parameter of the Gaussian.
      real(dp), intent(in) :: b

      !> Fine grid index each output bin maps to.
      integer, intent(out) :: gridIndex(0:M - 1)

      !> Reciprocal of the window's Fourier transform at each output bin.
      real(dp), intent(out) :: factor(0:M - 1)

      integer :: a      ! output bin
      integer :: k      ! the signed mode it holds

      do a = 0, M - 1
         k = mode_of_bin(a, M)
         gridIndex(a) = modulo(k, n)
         factor(a) = exp(b*(pi*real(k, dp)/real(n, dp))**2)
      end do

   end subroutine deconvolution_table

   ! =====================================================================
   !  Adjoint transforms: particles to modes
   ! =====================================================================

   !> Adjoint transform in three dimensions: structure factors of N charges.
   subroutine adjoint_nfft_3d(nParticle, t1, t2, t3, q, M1, M2, M3, m, Shat, phaseTimes)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> First coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t1(nParticle)

      !> Second coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t2(nParticle)

      !> Third coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t3(nParticle)

      !> Charge carried by each node.
      real(dp), intent(in) :: q(nParticle)

      !> Number of retained modes along the first axis.
      integer, intent(in) :: M1

      !> Number of retained modes along the second axis.
      integer, intent(in) :: M2

      !> Number of retained modes along the third axis.
      integer, intent(in) :: M3

      !> Stencil half-width.
      integer, intent(in) :: m

      !> Structure factors, in FFT bin order along all three axes.
      complex(dp), intent(out) :: Shat(0:M1 - 1, 0:M2 - 1, 0:M3 - 1)

      !> Wall time in seconds of the three steps, as
      !> [spread, transform, deconvolve].  For profiling only.
      real(dp), intent(out), optional :: phaseTimes(3)

      integer  :: n1, n2, n3
      integer  :: j
      integer  :: o1, o2, o3
      integer  :: c1, c2, c3
      integer  :: i1, i2, i3
      integer  :: a1, a2, a3
      real(dp) :: b
      real(dp) :: chargeTimesWeight
      real(dp) :: w1(-m:m), w2(-m:m), w3(-m:m)
      complex(dp), allocatable :: fineGrid(:, :, :)
      integer, allocatable  :: index1(:), index2(:), index3(:)
      real(dp), allocatable :: factor1(:), factor2(:), factor3(:)
      integer(8) :: tick0, tick1, tick2, tick3, tickRate

      n1 = oversampling*M1
      n2 = oversampling*M2
      n3 = oversampling*M3
      b = window_shape(m)
      allocate (fineGrid(0:n1 - 1, 0:n2 - 1, 0:n3 - 1))
      fineGrid = (0.0_dp, 0.0_dp)

      ! Step 1: spread.  Shared grid, overlapping stencils, atomic updates on
      ! the real part only.
      call system_clock(tick0, tickRate)
      !$omp parallel do default(shared) &
      !$omp private(j, c1, c2, c3, w1, w2, w3, o1, o2, o3, i1, i2, i3, chargeTimesWeight) &
      !$omp schedule(guided)
      do j = 1, nParticle
         if (q(j) == 0.0_dp) cycle
         call stencil_weights(t1(j), n1, m, b, c1, w1)
         call stencil_weights(t2(j), n2, m, b, c2, w2)
         call stencil_weights(t3(j), n3, m, b, c3, w3)
         do o3 = -m, m
            i3 = modulo(c3 + o3, n3)
            do o2 = -m, m
               i2 = modulo(c2 + o2, n2)
               chargeTimesWeight = q(j)*w2(o2)*w3(o3)
               do o1 = -m, m
                  i1 = modulo(c1 + o1, n1)
                  !$omp atomic update
                  fineGrid(i1, i2, i3)%re = fineGrid(i1, i2, i3)%re &
                                            + chargeTimesWeight*w1(o1)
               end do
            end do
         end do
      end do
      !$omp end parallel do
      call system_clock(tick1)

      ! Step 2: one ordinary transform over all three axes.
      call fft_3d_inplace(fineGrid, n1, n2, n3, 1)
      call system_clock(tick2)

      ! Step 3: keep the low modes and divide out the window.
      allocate (index1(0:M1 - 1), factor1(0:M1 - 1), index2(0:M2 - 1), factor2(0:M2 - 1), &
                index3(0:M3 - 1), factor3(0:M3 - 1))
      call deconvolution_table(M1, n1, b, index1, factor1)
      call deconvolution_table(M2, n2, b, index2, factor2)
      call deconvolution_table(M3, n3, b, index3, factor3)
      do a3 = 0, M3 - 1
         do a2 = 0, M2 - 1
            do a1 = 0, M1 - 1
               Shat(a1, a2, a3) = fineGrid(index1(a1), index2(a2), index3(a3)) &
                                  *factor1(a1)*factor2(a2)*factor3(a3)
            end do
         end do
      end do
      call system_clock(tick3)

      if (present(phaseTimes)) &
         phaseTimes = [real(tick1 - tick0, dp), real(tick2 - tick1, dp), &
                       real(tick3 - tick2, dp)]/real(tickRate, dp)

   end subroutine adjoint_nfft_3d

   ! =====================================================================
   !  Forward transforms: modes to particles
   ! =====================================================================

   !> Forward transform in three dimensions, the exact transpose of
   !> adjoint_nfft_3d.
   subroutine forward_nfft_3d(nParticle, t1, t2, t3, M1, M2, M3, m, fhat, f)

      !> Number of nodes to evaluate at.
      integer, intent(in) :: nParticle

      !> First coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t1(nParticle)

      !> Second coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t2(nParticle)

      !> Third coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t3(nParticle)

      !> Number of modes along the first axis.
      integer, intent(in) :: M1

      !> Number of modes along the second axis.
      integer, intent(in) :: M2

      !> Number of modes along the third axis.
      integer, intent(in) :: M3

      !> Stencil half-width.
      integer, intent(in) :: m

      !> Fourier coefficients, in FFT bin order along all three axes.
      complex(dp), intent(in) :: fhat(0:M1 - 1, 0:M2 - 1, 0:M3 - 1)

      !> The polynomial evaluated at each node.
      complex(dp), intent(out) :: f(nParticle)

      integer  :: n1, n2, n3
      integer  :: j, o1, o2, o3, c1, c2, c3, i1, i2, i3, a1, a2, a3
      real(dp) :: b
      complex(dp) :: accumulator
      complex(dp) :: outerWeight          ! product of the two outer weights
      real(dp) :: w1(-m:m), w2(-m:m), w3(-m:m)
      complex(dp), allocatable :: fineGrid(:, :, :)
      integer, allocatable  :: index1(:), index2(:), index3(:)
      real(dp), allocatable :: factor1(:), factor2(:), factor3(:)

      n1 = oversampling*M1
      n2 = oversampling*M2
      n3 = oversampling*M3
      b = window_shape(m)
      allocate (fineGrid(0:n1 - 1, 0:n2 - 1, 0:n3 - 1))
      fineGrid = (0.0_dp, 0.0_dp)

      ! Step 1: embed the coefficients, pre-divided by the window's transform.
      allocate (index1(0:M1 - 1), factor1(0:M1 - 1), index2(0:M2 - 1), factor2(0:M2 - 1), &
                index3(0:M3 - 1), factor3(0:M3 - 1))
      call deconvolution_table(M1, n1, b, index1, factor1)
      call deconvolution_table(M2, n2, b, index2, factor2)
      call deconvolution_table(M3, n3, b, index3, factor3)
      do a3 = 0, M3 - 1
         do a2 = 0, M2 - 1
            do a1 = 0, M1 - 1
               fineGrid(index1(a1), index2(a2), index3(a3)) = fhat(a1, a2, a3) &
                                                              *factor1(a1)*factor2(a2)*factor3(a3)
            end do
         end do
      end do

      ! Step 2: transform with the sign conjugate to the adjoint's.
      call fft_3d_inplace(fineGrid, n1, n2, n3, -1)

      ! Step 3: gather back to the nodes; one output per iteration, no atomics.
      !$omp parallel do default(shared) &
      !$omp private(j, c1, c2, c3, w1, w2, w3, o1, o2, o3, i1, i2, i3) &
      !$omp private(accumulator, outerWeight) schedule(guided)
      do j = 1, nParticle
         call stencil_weights(t1(j), n1, m, b, c1, w1)
         call stencil_weights(t2(j), n2, m, b, c2, w2)
         call stencil_weights(t3(j), n3, m, b, c3, w3)
         accumulator = (0.0_dp, 0.0_dp)
         do o3 = -m, m
            i3 = modulo(c3 + o3, n3)
            do o2 = -m, m
               i2 = modulo(c2 + o2, n2)
               outerWeight = cmplx(w2(o2)*w3(o3), 0.0_dp, dp)
               do o1 = -m, m
                  i1 = modulo(c1 + o1, n1)
                  accumulator = accumulator + w1(o1)*outerWeight*fineGrid(i1, i2, i3)
               end do
            end do
         end do
         f(j) = accumulator
      end do
      !$omp end parallel do

   end subroutine forward_nfft_3d

   ! =====================================================================
   !  Energy-only spreading onto a real grid
   ! =====================================================================

   !> Spread charges onto a padded real grid in three dimensions.
   subroutine spread_charges_real_3d(nParticle, t1, t2, t3, q, n1, n2, n3, m, grid)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> First coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t1(nParticle)

      !> Second coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t2(nParticle)

      !> Third coordinate of each node, on the unit torus.
      real(dp), intent(in) :: t3(nParticle)

      !> Charge carried by each node.
      real(dp), intent(in) :: q(nParticle)

      !> Logical extent of the fine grid along the first axis.
      integer, intent(in) :: n1

      !> Extent of the fine grid along the second axis.
      integer, intent(in) :: n2

      !> Extent of the fine grid along the third axis.
      integer, intent(in) :: n3

      !> Stencil half-width.
      integer, intent(in) :: m

      !> Padded real grid, accumulated into.
      real(dp), intent(inout) :: grid(0:2*(n1/2 + 1) - 1, 0:n2 - 1, 0:n3 - 1)

      integer  :: j, o1, o2, o3, c1, c2, c3, i1, i2, i3
      real(dp) :: b
      real(dp) :: chargeTimesWeight
      real(dp) :: w1(-m:m), w2(-m:m), w3(-m:m)

      b = window_shape(m)

      !$omp parallel do default(shared) &
      !$omp private(j, c1, c2, c3, w1, w2, w3, o1, o2, o3, i1, i2, i3, chargeTimesWeight) &
      !$omp schedule(guided)
      do j = 1, nParticle
         if (q(j) == 0.0_dp) cycle
         call stencil_weights(t1(j), n1, m, b, c1, w1)
         call stencil_weights(t2(j), n2, m, b, c2, w2)
         call stencil_weights(t3(j), n3, m, b, c3, w3)
         do o3 = -m, m
            i3 = modulo(c3 + o3, n3)
            do o2 = -m, m
               i2 = modulo(c2 + o2, n2)
               chargeTimesWeight = q(j)*w2(o2)*w3(o3)
               do o1 = -m, m
                  i1 = modulo(c1 + o1, n1)
                  !$omp atomic update
                  grid(i1, i2, i3) = grid(i1, i2, i3) + chargeTimesWeight*w1(o1)
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine spread_charges_real_3d

end module nfft

module ewald_fft_3d_parameters
   !> Parameter selection for the fast three-dimensional method. The real-space
   !> cutoff is tied to the mean particle spacing. Error estimates then determine
   !> the splitting parameter, Fourier cutoff and interpolation window, and the
   !> cell geometry determines the transform sizes.
   use ewald_constants, only: dp, pi
   use ewald_geometry, only: TCell3d
   use ewald_truncation, only: splitting_from_real_budget, cutoff_from_fourier_budget
   use nfft, only: window_cutoff_from_budget, oversampling
   use fft_backend, only: next_fast_length, minTransformLength
   implicit none

   private
   public :: TEwaldParameters3d, choose_parameters, grid_memory_gigabytes

   !> Calibrated factor in r_cut = realCutoffPrefactor * (V/N)^(1/3).
   !> Increasing it moves work from the grid to the pair sum.
   real(dp), parameter :: realCutoffPrefactor = 3.0_dp

   !> Fraction of the requested accuracy assigned to real-space truncation.
   real(dp), parameter :: realCutoffBudgetFraction = 0.2_dp

   !> Everything the fast method needs to evaluate one configuration.
   type :: TEwaldParameters3d

      !> Ewald splitting parameter, in inverse length.
      real(dp) :: alpha = 0.0_dp

      !> Real-space cutoff radius.
      real(dp) :: r_cut = 0.0_dp

      !> Fourier cutoff, as a physical wavenumber.
      real(dp) :: k_cut = 0.0_dp

      !> Half-width of the non-uniform transform's window stencil, in grid
      !> points.
      integer :: windowCutoff = 0

      !> Number of retained modes along each of the three axes.
      integer :: nModes(3) = 0

   end type TEwaldParameters3d

contains

   !> Work the whole parameter chain through for one configuration.
   !>
   !> An override is respected downstream: a caller who pins the cutoff gets a
   !> splitting parameter derived from that cutoff, not the default one, which
   !> is what lets a calibration sweep vary one parameter while the rest stay
   !> consistent with it.
   subroutine choose_parameters(cell, chargeSquareSum, nParticle, tolerance, params, &
                                alpha_in, r_cut_in, k_cut_in, nModes_in, windowCutoff_in)

      !> Geometry of the cell.
      type(TCell3d), intent(in) :: cell

      !> Sum of the squared charges, sum_i q_i^2.
      real(dp), intent(in) :: chargeSquareSum

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Requested accuracy of the whole evaluation.
      real(dp), intent(in) :: tolerance

      !> The resulting parameter set.
      type(TEwaldParameters3d), intent(out) :: params

      !> Override for the splitting parameter.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the Fourier cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      !> Override for the mode count, applied to all three axes at once.
      integer, intent(in), optional :: nModes_in

      !> Override for the window stencil half-width.
      integer, intent(in), optional :: windowCutoff_in

      real(dp) :: meanSpacing      ! (V/N)^(1/3), the typical distance between charges
      real(dp) :: cellListLimit    ! widest cutoff the linked-cell scheme allows
      integer  :: iAxis

      meanSpacing = (cell%volume/real(nParticle, dp))**(1.0_dp/3.0_dp)

      ! --- 1. real-space cutoff --------------------------------------------
      if (present(r_cut_in)) then
         params%r_cut = r_cut_in
      else if (present(alpha_in)) then
         ! A pinned splitting parameter implies a cutoff: the one at which the
         ! screened pair term has fallen to the tolerance.
         params%r_cut = sqrt(-log(tolerance))/alpha_in
      else
         params%r_cut = realCutoffPrefactor*meanSpacing

         ! The linked-cell scheme needs at least three boxes along every axis,
         ! each at least as wide as the cutoff; a wider cutoff would silently
         ! drop the short range back to a quadratic image sum.  The splitting
         ! parameter is rebalanced against whatever cutoff survives the cap, and
         ! the energy does not depend on the splitting, so capping costs no
         ! accuracy.  It is applied only when it still leaves a cutoff well
         ! above the mean spacing: for a genuinely small cell the image sum is
         ! valid and inexpensive.
         cellListLimit = 0.999_dp/(3.0_dp*maxval(cell%recLengths))
         if (cellListLimit >= 1.5_dp*meanSpacing) &
            params%r_cut = min(params%r_cut, cellListLimit)
      end if

      ! --- 2. splitting parameter ------------------------------------------
      if (present(alpha_in)) then
         params%alpha = alpha_in
      else
         params%alpha = splitting_from_real_budget(params%r_cut, chargeSquareSum, &
                                                   nParticle, cell%volume, &
                                                   realCutoffBudgetFraction*tolerance)
      end if
      if (params%alpha <= 0.0_dp) error stop "ewald_fft_3d: alpha must be positive"

      ! --- 3. Fourier cutoff -------------------------------------------------
      ! The estimate is written for a cubic box, so the cell is represented by
      ! the cube of the same volume.  The energy is not sensitive to that
      ! choice: it only shifts where the mode set is truncated, and the
      ! per-axis mode counts below use the true cell shape.
      params%k_cut = cutoff_from_fourier_budget(chargeSquareSum, nParticle, &
                                                cell%volume**(1.0_dp/3.0_dp), &
                                                params%alpha, tolerance)
      if (present(k_cut_in)) params%k_cut = k_cut_in

      ! --- 4. window stencil ------------------------------------------------
      params%windowCutoff = window_cutoff_from_budget(tolerance)
      if (present(windowCutoff_in)) params%windowCutoff = windowCutoff_in

      ! --- 5. mode counts ---------------------------------------------------
      do iAxis = 1, 3
         params%nModes(iAxis) = modes_along_axis(params%k_cut, cell%latLengths(iAxis))
      end do
      if (present(nModes_in)) params%nModes = nModes_in

   end subroutine choose_parameters

   !> Number of modes to keep along one axis so that the mode set contains
   !> every wavevector inside the cutoff.
   !>
   !> A mode index along axis i is m_i = k . a_i / (2*pi), so the cutoff sphere
   !> reaches k_cut*|a_i|/(2*pi) along that index, set by the length of the
   !> lattice vector and not by the perpendicular width of the cell.  Using the
   !> width instead undersizes the mode set for a sheared cell.  The index runs
   !> over both signs, hence the factor of two.
   pure function modes_along_axis(k_cut, latticeLength) result(nModes)

      !> Fourier cutoff, as a physical wavenumber.
      real(dp), intent(in) :: k_cut

      !> Length |a_i| of the lattice vector along this axis.
      real(dp), intent(in) :: latticeLength

      !> Number of modes to retain along this axis.
      integer :: nModes

      integer :: highestIndex   ! largest mode index inside the cutoff, plus a margin

      highestIndex = ceiling(k_cut*latticeLength/(2.0_dp*pi)) + 1
      nModes = max(next_fast_length(2*highestIndex), minTransformLength)

   end function modes_along_axis

   !> Peak grid memory of one evaluation, in gigabytes, predicted through the
   !> same parameter chain the evaluation itself uses.  Only the grids are
   !> counted; the cell list and the particle arrays are linear in the particle
   !> number and negligible beside them.
   !>
   !> The energy path allocates a single padded real grid, which the real-input
   !> transform overwrites with its half spectrum.  The potential and force
   !> path cannot use that shortcut and holds the full complex fine grid
   !> together with two mode arrays, roughly two and a half times as much.
   function grid_memory_gigabytes(cell, chargeSquareSum, nParticle, tolerance, &
                                  potentialForcePath, alpha_in, r_cut_in, k_cut_in, &
                                  nModes_in, windowCutoff_in) result(gigabytes)

      !> Geometry of the cell.
      type(TCell3d), intent(in) :: cell

      !> Sum of the squared charges.
      real(dp), intent(in) :: chargeSquareSum

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Requested accuracy.
      real(dp), intent(in) :: tolerance

      !> Predict the potential and force path instead of the energy path.
      logical, intent(in), optional :: potentialForcePath

      !> Override for the splitting parameter.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the Fourier cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      !> Override for the mode count on all three axes.
      integer, intent(in), optional :: nModes_in

      !> Override for the window stencil half-width.
      integer, intent(in), optional :: windowCutoff_in

      !> Predicted peak grid memory, in gigabytes.
      real(dp) :: gigabytes

      type(TEwaldParameters3d) :: params
      integer  :: n1, n2, n3        ! fine grid extents
      real(dp) :: nElements         ! number of array entries at the peak
      logical  :: forcePath

      call choose_parameters(cell, chargeSquareSum, nParticle, tolerance, params, &
                             alpha_in, r_cut_in, k_cut_in, nModes_in, windowCutoff_in)

      n1 = oversampling*params%nModes(1)
      n2 = oversampling*params%nModes(2)
      n3 = oversampling*params%nModes(3)

      forcePath = .false.
      if (present(potentialForcePath)) forcePath = potentialForcePath

      if (forcePath) then
         ! Complex fine grid, plus the coefficient array and the temporary that
         ! coexist with it while the force components are transformed.
         nElements = real(n1, dp)*real(n2, dp)*real(n3, dp) &
                     + 2.0_dp*real(params%nModes(1), dp)*real(params%nModes(2), dp) &
                     *real(params%nModes(3), dp)
         gigabytes = nElements*16.0_dp/(1024.0_dp**3)
      else
         ! One padded real grid; the padding of the first axis is what makes
         ! the real-input transform in place.
         nElements = real(2*(n1/2 + 1), dp)*real(n2, dp)*real(n3, dp)
         gigabytes = nElements*8.0_dp/(1024.0_dp**3)
      end if

   end function grid_memory_gigabytes

end module ewald_fft_3d_parameters

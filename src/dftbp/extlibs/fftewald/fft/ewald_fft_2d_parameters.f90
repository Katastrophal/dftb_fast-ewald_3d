module ewald_fft_2d_parameters
   !> Parameter selection for the fast two-dimensional method. The in-plane
   !> cutoffs follow the particle density and error estimates. The extent along
   !> the open direction determines whether the coplanar path applies and, when
   !> needed, the regularisation order, period and third transform size.
   use ewald_constants, only: dp, pi
   use ewald_geometry, only: TCell2d
   use ewald_truncation, only: splitting_from_real_budget, cutoff_from_fourier_budget
   use nfft, only: window_cutoff_from_budget, oversampling
   use fft_backend, only: next_fast_length, minTransformLength
   implicit none

   private
   public :: TEwaldParameters2d, choose_parameters, grid_memory_gigabytes
   public :: measure_extent, measure_thickness, is_monolayer
   public :: boundaryLayerWidth, set_boundary_layer_width, monolayerThreshold

   !> Threshold used to identify an exactly coplanar configuration.
   real(dp), parameter :: monolayerThreshold = 1.0e-10_dp

   !> Calibrated safety factor applied to both truncation budgets.
   real(dp), parameter :: truncationSafety = 0.04_dp

   !> Calibrated factor in r_cut = prefactor * sqrt(A/N) for a monolayer.
   real(dp), parameter :: monolayerCutoffPrefactor = 3.0_dp

   !> Corresponding calibrated factor for configurations of finite thickness.
   real(dp), parameter :: slabCutoffPrefactor = 22.0_dp

   !> Calibration factor of the boundary-layer resolution requirement: how many
   !> samples the mode set must place across the layer per matched derivative.
   real(dp), parameter :: layerResolutionFactor = 0.8_dp

   !> Calibration factor of the coverage requirement along the open direction.
   !> One leaves the nominal coverage in place; raising it enlarges the mode set
   !> where the residual tail there has to be pushed further below the request.
   real(dp), parameter :: normalCoverageFactor = 1.0_dp

   !> Relative width of the regularisation's boundary layer, as a fraction of
   !> the imposed period.
   !>
   !> A module variable rather than a constant only so that the calibration
   !> sweep can vary it.  Both the mode count along the open direction and the
   !> imposed period depend on it, so a caller who changes it must change it
   !> back afterwards; everything else assumes the shipped value.
   real(dp), protected :: boundaryLayerWidth = 0.1_dp

   !> Everything the fast method needs to evaluate one configuration.
   type :: TEwaldParameters2d

      !> Ewald splitting parameter, in inverse length.
      real(dp) :: alpha = 0.0_dp

      !> Real-space cutoff radius.
      real(dp) :: r_cut = 0.0_dp

      !> In-plane Fourier cutoff, as a physical wavenumber.
      real(dp) :: k_cut = 0.0_dp

      !> Half-width of the non-uniform transform's window stencil.
      integer :: windowCutoff = 0

      !> Smoothness of the kernel regularisation: the number of derivatives
      !> matched at each end of the boundary layer, counting the value.
      integer :: smoothness = 0

      !> Number of retained modes: the two in-plane axes and the open one.  The
      !> third entry is one for a monolayer.
      integer :: nModes(3) = 0

      !> Extent of the configuration along the open direction.
      real(dp) :: thickness = 0.0_dp

      !> Midpoint of that extent, about which the imposed period is centred.
      real(dp) :: centreZ = 0.0_dp

      !> Period imposed along the open direction by the regularisation, chosen
      !> so that every separation that can occur falls in the region where the
      !> regularised kernel equals the true one.  Zero for a monolayer, which
      !> imposes no period.
      real(dp) :: period = 0.0_dp

      !> Relative width of the boundary layer this parameter set was built
      !> with, carried so that the Fourier branch regularises with the same
      !> value the mode count along the open direction was sized for.
      real(dp) :: layerWidth = 0.0_dp

      !> Whether the configuration is coplanar and takes the single-mode path.
      logical :: monolayer = .true.

   end type TEwaldParameters2d

contains

   !> Override the width of the regularisation's boundary layer.  For
   !> calibration only; production callers should leave the shipped value
   !> alone.
   subroutine set_boundary_layer_width(width)

      !> New relative width, strictly between zero and one half.
      real(dp), intent(in) :: width

      if (width <= 0.0_dp .or. width >= 0.5_dp) &
         error stop "ewald_fft_2d: the boundary layer width must lie in (0,1/2)"
      boundaryLayerWidth = width

   end subroutine set_boundary_layer_width

   !> Extent of the configuration along the open direction, and its midpoint.
   !>
   !> Both quantities are computed here and stored in
   !> the parameter set, so that the parameters, the branch taken and the
   !> transform can never disagree about the geometry.
   pure subroutine measure_extent(positions, nParticle, thickness, centreZ)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Distance between the outermost charges along the open direction.
      real(dp), intent(out) :: thickness

      !> Midpoint between them.
      real(dp), intent(out) :: centreZ

      real(dp) :: lowestZ, highestZ

      lowestZ = minval(positions(1:nParticle, 3))
      highestZ = maxval(positions(1:nParticle, 3))
      thickness = highestZ - lowestZ
      centreZ = 0.5_dp*(lowestZ + highestZ)

   end subroutine measure_extent

   !> Extent alone, for callers that do not need the midpoint.
   pure function measure_thickness(positions, nParticle) result(thickness)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Distance between the outermost charges along the open direction.
      real(dp) :: thickness

      real(dp) :: centreZ

      call measure_extent(positions, nParticle, thickness, centreZ)

   end function measure_thickness

   !> Whether a configuration of the given thickness counts as coplanar.
   pure function is_monolayer(thickness) result(coplanar)

      !> Extent along the open direction, from measure_thickness.
      real(dp), intent(in) :: thickness

      !> True if the single-mode path applies.
      logical :: coplanar

      coplanar = (thickness <= monolayerThreshold)

   end function is_monolayer

   !> Work the whole parameter chain through for one configuration.
   !>
   !> As in the fully periodic case, an override is respected downstream: a
   !> caller who pins the cutoff gets a splitting parameter derived from it.
   subroutine choose_parameters(cell, chargeSquareSum, nParticle, thickness, centreZ, &
                                tolerance, params, alpha_in, r_cut_in, k_cut_in, &
                                normalModes_in, smoothness_in, windowCutoff_in)

      !> Geometry of the in-plane cell.
      type(TCell2d), intent(in) :: cell

      !> Sum of the squared charges, sum_i q_i^2.
      real(dp), intent(in) :: chargeSquareSum

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Extent of the configuration along the open direction.
      real(dp), intent(in) :: thickness

      !> Midpoint of that extent, from measure_extent.
      real(dp), intent(in) :: centreZ

      !> Requested accuracy of the whole evaluation.
      real(dp), intent(in) :: tolerance

      !> The resulting parameter set.
      type(TEwaldParameters2d), intent(out) :: params

      !> Override for the splitting parameter.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the in-plane Fourier cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      !> Override for the mode count along the open direction.
      integer, intent(in), optional :: normalModes_in

      !> Override for the regularisation smoothness.
      integer, intent(in), optional :: smoothness_in

      !> Override for the window stencil half-width.
      integer, intent(in), optional :: windowCutoff_in

      real(dp) :: meanSpacing      ! sqrt(A/N), the typical in-plane distance
      real(dp) :: cellListLimit    ! widest cutoff the linked-cell scheme allows
      integer  :: highestIndex

      params%thickness = thickness
      params%centreZ = centreZ
      params%layerWidth = boundaryLayerWidth
      params%monolayer = is_monolayer(thickness)
      meanSpacing = sqrt(cell%area/real(nParticle, dp))

      ! --- 1. real-space cutoff ---------------------------------------------
      if (present(r_cut_in)) then
         params%r_cut = r_cut_in
      else if (present(alpha_in)) then
         params%r_cut = sqrt(-log(tolerance))/alpha_in
      else
         if (params%monolayer) then
            params%r_cut = monolayerCutoffPrefactor*meanSpacing
         else
            params%r_cut = slabCutoffPrefactor*meanSpacing
         end if

         ! Cap at the widest cutoff the linked-cell scheme can handle, as in the
         ! fully periodic case: without it an elongated in-plane cell can push
         ! the cutoff past a third of its narrow width and drop the short range
         ! back to an explicit image sum.  The splitting parameter is rebalanced
         ! against whatever cutoff survives.  The cap is skipped for a genuinely
         ! small cell, where the image sum is valid and inexpensive.
         cellListLimit = 0.999_dp/(3.0_dp*maxval(cell%recLengths))
         if (cellListLimit >= 1.5_dp*meanSpacing) &
            params%r_cut = min(params%r_cut, cellListLimit)
      end if

      ! --- 2. splitting parameter -------------------------------------------
      ! The cubic box equivalent to a two-dimensional cell is taken to have the
      ! in-plane cell's edge, so its volume is the area raised to three halves.
      ! That is well defined even for a monolayer, whose physical volume is
      ! zero, and the energy does not depend on the choice: it only shifts work
      ! between the two branches.
      if (present(alpha_in)) then
         params%alpha = alpha_in
      else
         params%alpha = splitting_from_real_budget(params%r_cut, chargeSquareSum, &
                                                   nParticle, cell%area**1.5_dp, &
                                                   truncationSafety*tolerance)
      end if
      if (params%alpha <= 0.0_dp) error stop "ewald_fft_2d: alpha must be positive"

      ! --- 3. in-plane Fourier cutoff ---------------------------------------
      params%k_cut = cutoff_from_fourier_budget(chargeSquareSum, nParticle, &
                                                sqrt(cell%area), params%alpha, &
                                                truncationSafety*tolerance)
      if (present(k_cut_in)) params%k_cut = k_cut_in

      ! --- 4. window stencil -------------------------------------------------
      params%windowCutoff = window_cutoff_from_budget(tolerance)
      if (present(windowCutoff_in)) params%windowCutoff = windowCutoff_in

      ! --- 5. regularisation smoothness --------------------------------------
      ! The Hermite closure has no error estimate of its own. Its observed
      ! error decay in the number of matched derivatives is at least as fast as
      ! the window's decay in the stencil width, so the same budget-derived
      ! value serves both.
      params%smoothness = window_cutoff_from_budget(tolerance)
      if (present(smoothness_in)) params%smoothness = smoothness_in

      ! --- 6a. in-plane mode counts -------------------------------------------
      ! As in three dimensions, an index reaches k_cut*|a_i|/(2*pi) along axis
      ! i, set by the lattice vector length and not by the perpendicular width
      ! of the cell.  The two agree only for an orthogonal cell; using the
      ! width instead undersizes the mode set for an oblique one.
      highestIndex = ceiling(params%k_cut*cell%latLengths(1)/(2.0_dp*pi)) + 1
      params%nModes(1) = max(next_fast_length(2*highestIndex), minTransformLength)
      highestIndex = ceiling(params%k_cut*cell%latLengths(2)/(2.0_dp*pi)) + 1
      params%nModes(2) = max(next_fast_length(2*highestIndex), minTransformLength)

      ! --- 6b. modes along the open direction ---------------------------------
      ! Two requirements, and the mode count has to meet both.  Resolution: the
      ! boundary layer occupies a fixed fraction of the period and the
      ! polynomial living in it has degree 2p-1, so the mode set must place
      ! enough samples across the layer to represent it.  That scales with the
      ! smoothness and not with the thickness, which is why it is a floor.
      ! Coverage: the same physical cutoff has to be covered along the open
      ! direction as in the plane, a mode index l mapping to 2*pi*l/period.
      !
      ! A monolayer needs neither: it has a single mode there.
      if (params%monolayer) then
         params%period = 0.0_dp
         params%nModes(3) = 1
      else
         params%period = thickness/(0.5_dp - boundaryLayerWidth)
         params%nModes(3) = next_fast_length( &
                            max(ceiling(layerResolutionFactor*real(params%smoothness, dp) &
                                        /boundaryLayerWidth), &
                                ceiling(normalCoverageFactor*params%period*params%k_cut/pi)))
      end if
      if (present(normalModes_in)) params%nModes(3) = normalModes_in

   end subroutine choose_parameters

   !> Peak grid memory of one evaluation, in gigabytes, predicted through the
   !> same parameter chain the evaluation itself uses.  Which grid is predicted
   !> follows the computed thickness, by the criterion that selects the branch
   !> of the evaluation.
   function grid_memory_gigabytes(cell, chargeSquareSum, nParticle, thickness, centreZ, &
                                  tolerance, potentialForcePath, alpha_in, r_cut_in, &
                                  k_cut_in, normalModes_in, smoothness_in, &
                                  windowCutoff_in) result(gigabytes)

      !> Geometry of the in-plane cell.
      type(TCell2d), intent(in) :: cell

      !> Sum of the squared charges.
      real(dp), intent(in) :: chargeSquareSum

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Extent of the configuration along the open direction.
      real(dp), intent(in) :: thickness

      !> Midpoint of that extent, from measure_extent.
      real(dp), intent(in) :: centreZ

      !> Requested accuracy.
      real(dp), intent(in) :: tolerance

      !> Predict the potential and force path instead of the energy path.
      logical, intent(in), optional :: potentialForcePath

      !> Override for the splitting parameter.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the in-plane Fourier cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      !> Override for the mode count along the open direction.
      integer, intent(in), optional :: normalModes_in

      !> Override for the regularisation smoothness.
      integer, intent(in), optional :: smoothness_in

      !> Override for the window stencil half-width.
      integer, intent(in), optional :: windowCutoff_in

      !> Predicted peak grid memory, in gigabytes.
      real(dp) :: gigabytes

      type(TEwaldParameters2d) :: params
      integer  :: n1, n2, n3
      real(dp) :: modeCount        ! product of the three mode counts
      logical  :: forcePath

      call choose_parameters(cell, chargeSquareSum, nParticle, thickness, centreZ, &
                             tolerance, params, alpha_in, r_cut_in, k_cut_in, &
                             normalModes_in, smoothness_in, windowCutoff_in)

      n1 = oversampling*params%nModes(1)
      n2 = oversampling*params%nModes(2)
      ! A monolayer's fine grid has a single plane, not an oversampled one.
      n3 = 1
      if (.not. params%monolayer) n3 = oversampling*params%nModes(3)

      forcePath = .false.
      if (present(potentialForcePath)) forcePath = potentialForcePath

      modeCount = real(params%nModes(1), dp)*real(params%nModes(2), dp)

      if (forcePath) then
         ! The potential and force path works on signed structure factors and
         ! cannot use the half spectrum.  At its peak the full complex fine
         ! grid coexists with two complex mode arrays, all at sixteen bytes per
         ! entry.
         if (params%monolayer) then
            gigabytes = real(2 + oversampling*oversampling, dp)*modeCount &
                        *16.0_dp/(1024.0_dp**3)
         else
            gigabytes = real(2 + oversampling**3, dp)*modeCount &
                        *real(params%nModes(3), dp)*16.0_dp/(1024.0_dp**3)
         end if
      else
         ! One padded real grid, whose first dimension carries the padding the
         ! in-place real-input transform requires.
         gigabytes = real(2*(n1/2 + 1), dp)*real(n2, dp)*real(n3, dp) &
                     *8.0_dp/(1024.0_dp**3)
      end if

   end function grid_memory_gigabytes

end module ewald_fft_2d_parameters

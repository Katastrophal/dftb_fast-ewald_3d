module ewald_fft_2d
   !> Fast Ewald summation for a cell that is periodic in two directions and
   !> open in the third.
   !>
   !> The real-space branch uses a linked-cell decomposition. The Fourier branch
   !> uses a two-dimensional transform for coplanar configurations and a
   !> regularised three-dimensional kernel otherwise. This facade validates the
   !> arguments, selects the path and assembles the contributions. Only
   !> charge-neutral configurations are accepted.
   use ewald_constants, only: dp
   use ewald_geometry, only: TCell2d, cell_metrics_2d
   use ewald_validation, only: check_configuration
   use ewald_self, only: self_energy, self_potential
   use ewald_fft_2d_parameters, only: TEwaldParameters2d, choose_parameters, &
                                      grid_memory_gigabytes, measure_extent, &
                                      measure_thickness, is_monolayer, &
                                      monolayerThreshold, boundaryLayerWidth, &
                                      set_boundary_layer_width
   use ewald_fft_2d_real, only: short_range_energy, short_range_potential_force
   use ewald_fft_2d_fourier, only: long_range_energy, &
                                   long_range_potential_force_monolayer, &
                                   long_range_potential_force_slab
   implicit none

   private
   public :: ewald_energy, ewald_potential_force, ewald_parameters, ewald_grid_memory
   public :: TEwaldParameters2d
   !> Re-exported so that a caller deciding whether this module suits a
   !> geometry tests the same criterion the module branches on.
   public :: monolayerThreshold, measure_extent, measure_thickness, is_monolayer
   !> Re-exported for the calibration of the regularisation; read-only outside
   !> the setter.
   public :: boundaryLayerWidth, set_boundary_layer_width

   !> Accuracy used when the caller does not request one.
   real(dp), parameter :: defaultTolerance = 1.0e-10_dp

contains

   !> Total electrostatic energy of the two-dimensionally periodic cell.
   function ewald_energy(positions, charges, nParticle, latVecs, tol, &
                         alpha_in, r_cut_in, k_cut_in, normalModes_in, smoothness_in, &
                         windowCutoff_in, phaseTimes) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.  The third
      !> component is the coordinate along the open direction.
      real(dp), intent(in) :: positions(:, :)

      !> Charges, in units where the Coulomb prefactor is one.
      real(dp), intent(in) :: charges(:)

      !> Number of charges in the cell.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.  Rows 1 and 2 span the periodic plane;
      !> row 3 is ignored.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy.  Defaults to 1e-10.
      real(dp), intent(in), optional :: tol

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

      !> Wall time in seconds of the five phases of one evaluation, as
      !> [short range, spread, transform, extract, mode sum].  For profiling
      !> only.
      real(dp), intent(out), optional :: phaseTimes(5)

      !> Total energy of the cell.
      real(dp) :: energy

      type(TCell2d) :: cell
      type(TEwaldParameters2d) :: params
      real(dp) :: tolerance
      real(dp) :: chargeSquareSum
      real(dp) :: thickness, centreZ
      real(dp) :: shortRangeEnergy
      real(dp) :: longRangeEnergy
      real(dp) :: shortRangeTime
      real(dp) :: longRangeTimes(4)
      integer(8) :: tick0, tick1, tickRate

      call check_configuration("ewald_fft_2d", positions, charges, nParticle)

      cell = cell_metrics_2d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_fft_2d: tol must lie in (0,1)"

      chargeSquareSum = sum(charges(1:nParticle)**2)
      call measure_extent(positions, nParticle, thickness, centreZ)

      call choose_parameters(cell, chargeSquareSum, nParticle, thickness, centreZ, &
                             tolerance, params, alpha_in, r_cut_in, k_cut_in, &
                             normalModes_in, smoothness_in, windowCutoff_in)

      call system_clock(tick0, tickRate)
      shortRangeEnergy = short_range_energy(positions, charges, nParticle, cell, &
                                            params%alpha, params%r_cut)
      call system_clock(tick1)
      shortRangeTime = real(tick1 - tick0, dp)/real(tickRate, dp)

      longRangeEnergy = long_range_energy(positions, charges, nParticle, cell%recVecs, &
                                          cell%area, params, longRangeTimes)

      ! Three terms, not four: the zero in-plane mode is inside the Fourier
      ! branch rather than beside it, because the kernel carries it.
      energy = shortRangeEnergy + longRangeEnergy &
               + self_energy(charges, nParticle, params%alpha)

      if (present(phaseTimes)) phaseTimes = [shortRangeTime, longRangeTimes]

   end function ewald_energy

   !> Electrostatic potential at each charge and the force acting on it.
   !>
   !> The branch is chosen on the same criterion as the energy: a coplanar
   !> configuration takes the monolayer transform, anything thicker the slab
   !> one.  They differ in what the Fourier branch does, not in the splitting
   !> or the parameter chain.
   !>
   !> On a monolayer the force normal to the layer vanishes exactly: a coplanar
   !> charge distribution is symmetric under reflection in its own plane, so the
   !> potential is even in the normal coordinate.  A slab carries no such
   !> symmetry and its normal force is computed like the other two.
   !>
   !> The slab path has a caveat the energy does not.  Its kernel is smooth only
   !> to a finite order across the joint of the regularisation, so
   !> differentiating it for the force costs one order of that smoothness, and
   !> the shipped smoothness and mode-count floor are calibrated on the energy.
   !> The test suite checks the resulting force accuracy.
   !>
   !> The force argument is optional; omitting it skips the extra transforms.
   subroutine ewald_potential_force(positions, charges, nParticle, latVecs, tol, &
                                    pot, force, energy, alpha_in, r_cut_in, k_cut_in, &
                                    normalModes_in, smoothness_in, windowCutoff_in)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy.
      real(dp), intent(in), optional :: tol

      !> Potential at each charge, pot(i) = phi(x_i).  Size at least nParticle.
      real(dp), intent(out) :: pot(:)

      !> Force on each charge, force(:, i) = -dU/dx_i, shape (3, nParticle).
      real(dp), intent(out), optional :: force(:, :)

      !> Total energy, identical to what ewald_energy returns.
      real(dp), intent(out), optional :: energy

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

      type(TCell2d) :: cell
      type(TEwaldParameters2d) :: params
      real(dp) :: tolerance
      real(dp) :: chargeSquareSum
      real(dp) :: thickness, centreZ
      real(dp) :: shortRangeEnergy
      real(dp) :: longRangeEnergy
      real(dp), allocatable :: shortRangePot(:), longRangePot(:)
      real(dp), allocatable :: shortRangeForce(:, :), longRangeForce(:, :)
      integer :: i

      call check_configuration("ewald_fft_2d", positions, charges, nParticle)
      if (size(pot) < nParticle) error stop "ewald_fft_2d: pot must be sized (N)"
      if (present(force)) then
         if (size(force, 1) /= 3 .or. size(force, 2) < nParticle) &
            error stop "ewald_fft_2d: force must be sized (3,N)"
      end if

      cell = cell_metrics_2d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_fft_2d: tol must lie in (0,1)"

      chargeSquareSum = sum(charges(1:nParticle)**2)
      call measure_extent(positions, nParticle, thickness, centreZ)

      call choose_parameters(cell, chargeSquareSum, nParticle, thickness, centreZ, &
                             tolerance, params, alpha_in, r_cut_in, k_cut_in, &
                             normalModes_in, smoothness_in, windowCutoff_in)

      allocate (shortRangePot(nParticle), longRangePot(nParticle))
      allocate (shortRangeForce(3, nParticle), longRangeForce(3, nParticle))

      call short_range_potential_force(positions, charges, nParticle, cell, &
                                       params%alpha, params%r_cut, &
                                       shortRangePot, shortRangeForce)

      if (params%monolayer) then
         if (present(force)) then
            call long_range_potential_force_monolayer(positions, charges, nParticle, &
                                                      cell%recVecs, cell%area, params, &
                                                      longRangeEnergy, longRangePot, &
                                                      longRangeForce)
         else
            call long_range_potential_force_monolayer(positions, charges, nParticle, &
                                                      cell%recVecs, cell%area, params, &
                                                      longRangeEnergy, longRangePot)
         end if
      else
         if (present(force)) then
            call long_range_potential_force_slab(positions, charges, nParticle, &
                                                 cell%recVecs, cell%area, params, &
                                                 longRangeEnergy, longRangePot, &
                                                 longRangeForce)
         else
            call long_range_potential_force_slab(positions, charges, nParticle, &
                                                 cell%recVecs, cell%area, params, &
                                                 longRangeEnergy, longRangePot)
         end if
      end if

      ! Self term only.  It is the same at every atom and independent of the
      ! positions, so it carries no force, and there is nothing else to add:
      ! the uniform-sheet contribution sits inside the Fourier branch, in the
      ! zero in-plane mode.
      do i = 1, nParticle
         pot(i) = shortRangePot(i) + longRangePot(i) &
                  + self_potential(charges(i), params%alpha)
      end do

      if (present(force)) then
         do i = 1, nParticle
            force(:, i) = shortRangeForce(:, i) + longRangeForce(:, i)
         end do
      end if

      if (present(energy)) then
         shortRangeEnergy = 0.5_dp*dot_product(charges(1:nParticle), &
                                               shortRangePot(1:nParticle))
         energy = shortRangeEnergy + longRangeEnergy &
                  + self_energy(charges, nParticle, params%alpha)
      end if

   end subroutine ewald_potential_force

   !> The parameters this module would use for a given configuration, including
   !> the computed thickness.
   !>
   !> A calibration experiment has to pin one parameter and sweep another, which
   !> means it needs the shipped values of the rest.  There is one
   !> implementation of the chain and this exposes it.
   subroutine ewald_parameters(positions, charges, nParticle, latVecs, tol, params, &
                               alpha_in, r_cut_in, k_cut_in, normalModes_in, &
                               smoothness_in, windowCutoff_in)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy.
      real(dp), intent(in), optional :: tol

      !> The parameters an evaluation with these arguments would use.
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

      type(TCell2d) :: cell
      real(dp) :: tolerance
      real(dp) :: thickness, centreZ

      call check_configuration("ewald_fft_2d", positions, charges, nParticle)
      cell = cell_metrics_2d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol

      call measure_extent(positions, nParticle, thickness, centreZ)
      call choose_parameters(cell, sum(charges(1:nParticle)**2), nParticle, thickness, &
                             centreZ, tolerance, params, alpha_in, r_cut_in, k_cut_in, &
                             normalModes_in, smoothness_in, windowCutoff_in)

   end subroutine ewald_parameters

   !> Peak grid memory of one evaluation, in gigabytes.
   !>
   !> Only the coordinate along the open direction of the positions is read, to
   !> measure the thickness.  Omitting them describes a monolayer; a slab caller
   !> that omits them would be told the in-plane grid alone and would miss the
   !> mode count along the open direction, which is the larger factor by far.
   function ewald_grid_memory(charges, nParticle, latVecs, tol, positions, &
                              potentialForcePath, alpha_in, r_cut_in, k_cut_in, &
                              normalModes_in, smoothness_in, windowCutoff_in) &
      result(gigabytes)

      !> Charges.  Only their squared sum is read.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy.
      real(dp), intent(in), optional :: tol

      !> Positions, read only to measure the thickness.  Omitting them
      !> describes a monolayer.
      real(dp), intent(in), optional :: positions(:, :)

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

      type(TCell2d) :: cell
      real(dp) :: tolerance
      real(dp) :: thickness, centreZ

      cell = cell_metrics_2d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol

      thickness = 0.0_dp
      centreZ = 0.0_dp
      if (present(positions)) call measure_extent(positions, nParticle, thickness, centreZ)

      gigabytes = grid_memory_gigabytes(cell, sum(charges(1:nParticle)**2), nParticle, &
                                        thickness, centreZ, tolerance, &
                                        potentialForcePath, alpha_in, r_cut_in, &
                                        k_cut_in, normalModes_in, smoothness_in, &
                                        windowCutoff_in)

   end function ewald_grid_memory

end module ewald_fft_2d

module ewald_fft_3d
   !> Fast Ewald summation for a cell that is periodic in all three directions.
   !>
   !> The real-space branch uses a linked-cell decomposition, and the Fourier
   !> branch obtains the structure factors through a non-uniform transform. This
   !> facade validates the arguments, selects the parameters and assembles the
   !> contributions. Only charge-neutral configurations are accepted.
   use ewald_constants, only: dp
   use ewald_geometry, only: TCell3d, cell_metrics_3d
   use ewald_validation, only: check_configuration
   use ewald_self, only: self_energy, self_potential
   use ewald_fft_3d_parameters, only: TEwaldParameters3d, choose_parameters, &
                                      grid_memory_gigabytes
   use ewald_fft_3d_real, only: short_range_energy, short_range_potential_force
   use ewald_fft_3d_fourier, only: long_range_energy, long_range_potential_force
   implicit none

   private
   public :: ewald_energy, ewald_potential_force, ewald_parameters, ewald_grid_memory
   public :: TEwaldParameters3d

   !> Accuracy used when the caller does not request one.
   real(dp), parameter :: defaultTolerance = 1.0e-10_dp

contains

   !> Total electrostatic energy of the periodic cell.
   function ewald_energy(positions, charges, nParticle, latVecs, tol, &
                         alpha_in, r_cut_in, k_cut_in, nModes_in, windowCutoff_in, &
                         phaseTimes) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges, in units where the Coulomb prefactor is one.
      real(dp), intent(in) :: charges(:)

      !> Number of charges in the cell.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy.  Defaults to 1e-10.
      real(dp), intent(in), optional :: tol

      !> Override for the splitting parameter.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the Fourier cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      !> Override for the mode count, applied to all three axes.
      integer, intent(in), optional :: nModes_in

      !> Override for the window stencil half-width.
      integer, intent(in), optional :: windowCutoff_in

      !> Wall time in seconds of the five phases of one evaluation, as
      !> [short range, spread, transform, extract, mode sum].  For profiling
      !> only; the evaluation is unaffected by asking for it.
      real(dp), intent(out), optional :: phaseTimes(5)

      !> Total energy of the cell.
      real(dp) :: energy

      type(TCell3d) :: cell
      type(TEwaldParameters3d) :: params
      real(dp) :: tolerance
      real(dp) :: chargeSquareSum      ! sum_i q_i^2
      real(dp) :: shortRangeEnergy
      real(dp) :: longRangeEnergy
      real(dp) :: shortRangeTime
      real(dp) :: longRangeTimes(4)
      integer(8) :: tick0, tick1, tickRate

      call check_configuration("ewald_fft_3d", positions, charges, nParticle)

      cell = cell_metrics_3d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_fft_3d: tol must lie in (0,1)"

      chargeSquareSum = sum(charges(1:nParticle)**2)
      call choose_parameters(cell, chargeSquareSum, nParticle, tolerance, params, &
                             alpha_in, r_cut_in, k_cut_in, nModes_in, windowCutoff_in)

      call system_clock(tick0, tickRate)
      shortRangeEnergy = short_range_energy(positions, charges, nParticle, cell, &
                                            params%alpha, params%r_cut)
      call system_clock(tick1)
      shortRangeTime = real(tick1 - tick0, dp)/real(tickRate, dp)

      longRangeEnergy = long_range_energy(positions, charges, nParticle, cell%recVecs, &
                                          cell%volume, params%alpha, params%nModes, &
                                          params%windowCutoff, longRangeTimes)

      energy = shortRangeEnergy + longRangeEnergy &
               + self_energy(charges, nParticle, params%alpha)

      if (present(phaseTimes)) phaseTimes = [shortRangeTime, longRangeTimes]

   end function ewald_energy

   !> Electrostatic potential at each charge and the force acting on it: the
   !> quantities a self-consistent-charge electronic-structure code consumes.
   !> The splitting and the parameter chain are those of ewald_energy; what is
   !> added is the interpolation of the Fourier field back to the particles.
   !>
   !> The force argument is optional, and omitting it saves three of the four
   !> transforms, which is worth doing inside a self-consistency cycle that
   !> needs the potential at every iteration but the force only once.
   subroutine ewald_potential_force(positions, charges, nParticle, latVecs, tol, &
                                    pot, force, energy, alpha_in, r_cut_in, k_cut_in, &
                                    nModes_in, windowCutoff_in)

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
      !> This is the physical force; a code that stores energy gradients wants
      !> its negative.
      real(dp), intent(out), optional :: force(:, :)

      !> Total energy, identical to what ewald_energy returns.
      real(dp), intent(out), optional :: energy

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

      type(TCell3d) :: cell
      type(TEwaldParameters3d) :: params
      real(dp) :: tolerance
      real(dp) :: chargeSquareSum
      real(dp) :: shortRangeEnergy
      real(dp) :: longRangeEnergy
      real(dp), allocatable :: shortRangePot(:), longRangePot(:)
      real(dp), allocatable :: shortRangeForce(:, :), longRangeForce(:, :)
      integer :: i

      call check_configuration("ewald_fft_3d", positions, charges, nParticle)
      if (size(pot) < nParticle) error stop "ewald_fft_3d: pot must be sized (N)"
      if (present(force)) then
         if (size(force, 1) /= 3 .or. size(force, 2) < nParticle) &
            error stop "ewald_fft_3d: force must be sized (3,N)"
      end if

      cell = cell_metrics_3d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_fft_3d: tol must lie in (0,1)"

      chargeSquareSum = sum(charges(1:nParticle)**2)
      call choose_parameters(cell, chargeSquareSum, nParticle, tolerance, params, &
                             alpha_in, r_cut_in, k_cut_in, nModes_in, windowCutoff_in)

      allocate (shortRangePot(nParticle), longRangePot(nParticle))
      allocate (shortRangeForce(3, nParticle), longRangeForce(3, nParticle))

      call short_range_potential_force(positions, charges, nParticle, cell, &
                                       params%alpha, params%r_cut, &
                                       shortRangePot, shortRangeForce)

      if (present(force)) then
         call long_range_potential_force(positions, charges, nParticle, cell%recVecs, &
                                         cell%volume, params%alpha, params%nModes, &
                                         params%windowCutoff, longRangeEnergy, &
                                         longRangePot, longRangeForce)
      else
         call long_range_potential_force(positions, charges, nParticle, cell%recVecs, &
                                         cell%volume, params%alpha, params%nModes, &
                                         params%windowCutoff, longRangeEnergy, &
                                         longRangePot)
      end if

      ! The self term is the same at every atom and independent of the
      ! positions, so it contributes no force.
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
         ! The short-range energy is recovered exactly from its own potential,
         ! which is cheaper than summing the pairs a second time.
         shortRangeEnergy = 0.5_dp*dot_product(charges(1:nParticle), &
                                               shortRangePot(1:nParticle))
         energy = shortRangeEnergy + longRangeEnergy &
                  + self_energy(charges, nParticle, params%alpha)
      end if

   end subroutine ewald_potential_force

   !> The parameters this module would use for a given configuration.
   !>
   !> A calibration experiment has to pin one parameter and sweep another, which
   !> means it needs the shipped values of the rest.  There is one
   !> implementation of the chain and this exposes it, so that no driver has to
   !> re-derive it.
   subroutine ewald_parameters(positions, charges, nParticle, latVecs, tol, params, &
                               alpha_in, r_cut_in, k_cut_in, nModes_in, windowCutoff_in)

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
      type(TEwaldParameters3d), intent(out) :: params

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

      type(TCell3d) :: cell
      real(dp) :: tolerance

      call check_configuration("ewald_fft_3d", positions, charges, nParticle)
      cell = cell_metrics_3d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol

      call choose_parameters(cell, sum(charges(1:nParticle)**2), nParticle, tolerance, &
                             params, alpha_in, r_cut_in, k_cut_in, nModes_in, &
                             windowCutoff_in)

   end subroutine ewald_parameters

   !> Peak grid memory of one evaluation, in gigabytes, so that a caller
   !> sweeping the system size or the accuracy can gate on it before trying an
   !> evaluation that would not fit.  The prediction runs the same parameter
   !> chain the evaluation does.
   function ewald_grid_memory(charges, nParticle, latVecs, tol, potentialForcePath, &
                              alpha_in, r_cut_in, k_cut_in, nModes_in, windowCutoff_in) &
      result(gigabytes)

      !> Charges.  Only their squared sum is read, the positions playing no
      !> part in how large the grid is.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy.
      real(dp), intent(in), optional :: tol

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

      type(TCell3d) :: cell
      real(dp) :: tolerance

      cell = cell_metrics_3d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol

      gigabytes = grid_memory_gigabytes(cell, sum(charges(1:nParticle)**2), nParticle, &
                                        tolerance, potentialForcePath, alpha_in, &
                                        r_cut_in, k_cut_in, nModes_in, windowCutoff_in)

   end function ewald_grid_memory

end module ewald_fft_3d

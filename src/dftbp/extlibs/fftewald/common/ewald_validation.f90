module ewald_validation
   !> Input validation and numerical guards shared by all four methods.
   use ewald_constants, only: dp
   implicit none

   private
   public :: check_configuration
   public :: minSeparationSquared, minWavenumberSquared

   !> Squared separation below which two charges count as coincident.  A pair at
   !> zero distance has an infinite bare interaction, so this is normally an
   !> input error rather than something to be regularised away; the fast
   !> two-dimensional real-space branch is the one place that uses it as a
   !> filter instead, and says so at the use site.
   real(dp), parameter :: minSeparationSquared = 1.0e-24_dp

   !> Squared wavenumber below which a Fourier mode counts as the zero mode.
   real(dp), parameter :: minWavenumberSquared = 1.0e-24_dp

   !> Largest net charge, as a fraction of the total absolute charge, that still
   !> counts as neutral.
   real(dp), parameter :: neutralityTolerance = 1.0e-8_dp

contains

   !> Reject the argument shapes that would otherwise be read out of bounds or
   !> silently misinterpreted, and the cells no method here is defined for.
   subroutine check_configuration(caller, positions, charges, nParticle)
      character(len=*), intent(in) :: caller
      real(dp), intent(in) :: positions(:, :)
      real(dp), intent(in) :: charges(:)
      integer, intent(in) :: nParticle

      real(dp) :: netCharge, absoluteCharge

      if (nParticle <= 0) error stop caller//": nParticle must be positive"
      if (size(positions, 2) /= 3) error stop caller//": positions must be (N,3)"
      if (size(positions, 1) < nParticle) error stop caller//": too few positions"
      if (size(charges) < nParticle) error stop caller//": too few charges"

      ! Only a neutral cell has a finite energy without a compensating
      ! background, and no background is implemented.  The imbalance is judged
      ! against the total absolute charge rather than against a fixed number, so
      ! that the test means the same thing whatever the charges are scaled in.
      netCharge = sum(charges(1:nParticle))
      absoluteCharge = sum(abs(charges(1:nParticle)))
      if (abs(netCharge) > neutralityTolerance*absoluteCharge) &
         error stop caller//": the cell must be charge neutral"

   end subroutine check_configuration

end module ewald_validation

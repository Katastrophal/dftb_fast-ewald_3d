module ewald_self
   !> Analytic self-interaction corrections shared by all four methods.
   use ewald_constants, only: dp, sqrt_pi
   implicit none

   private
   public :: self_energy, self_potential

contains

   !> Self-interaction correction to the total energy.
   pure function self_energy(charges, nParticle, alpha) result(energy)
      real(dp), intent(in) :: charges(:)
      integer, intent(in) :: nParticle
      real(dp), intent(in) :: alpha
      real(dp) :: energy

      ! Slice, in case the caller passes an oversized array.
      energy = -alpha/sqrt_pi*sum(charges(1:nParticle)**2)

   end function self_energy

   !> Self-interaction correction to the potential at one charge.  It is twice
   !> the energy's share per charge, the energy carrying a factor of one half,
   !> so that U = (1/2) sum_i q_i phi(x_i) still holds.
   pure function self_potential(charge, alpha) result(potential)
      real(dp), intent(in) :: charge
      real(dp), intent(in) :: alpha
      real(dp) :: potential

      potential = -2.0_dp*alpha/sqrt_pi*charge

   end function self_potential

end module ewald_self

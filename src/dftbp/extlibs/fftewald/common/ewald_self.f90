module ewald_self
   !> Analytic self-interaction corrections shared by all four methods.
   use ewald_constants, only: dp, sqrt_pi
   implicit none

   private
   public :: self_potential

contains

   !> Self-interaction correction to the potential at one charge.
   pure function self_potential(charge, alpha) result(potential)
      real(dp), intent(in) :: charge
      real(dp), intent(in) :: alpha
      real(dp) :: potential

      potential = -2.0_dp*alpha/sqrt_pi*charge

   end function self_potential

end module ewald_self

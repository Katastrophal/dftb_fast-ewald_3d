module ewald_constants
   !> Working precision and the mathematical constants used throughout the
   !> library.
   use iso_fortran_env, only: real64
   implicit none

   private
   public :: dp, pi, sqrt_pi

   !> Working precision of every real quantity in the library.
   integer, parameter :: dp = real64

   !> Computed rather than typed, so that they carry the full precision of dp.
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: sqrt_pi = sqrt(pi)

end module ewald_constants

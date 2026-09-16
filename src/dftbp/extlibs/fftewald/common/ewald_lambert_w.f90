module ewald_lambert_w
   !> Principal branch of the Lambert W function, W(x) e^{W(x)} = x.
   !>
   !> The parameter recipes of both fast methods invert an error estimate of
   !> the form "constant times exp(-y) times a power of y equals the requested
   !> accuracy" for y, whose solution is a Lambert W value.  Only the principal
   !> branch and only non-negative arguments occur, the inverted quantity being
   !> an error budget.
   use ewald_constants, only: dp
   implicit none

   private
   public :: lambert_w0

   !> Relative size of the correction at which the iteration stops.
   real(dp), parameter :: converged = 1.0e-15_dp

   !> Safety net only; the loop exits on `converged` for every argument that
   !> occurs here.
   integer, parameter :: max_iterations = 60

contains

   !> Principal branch W0(x) for x >= 0, by Halley iteration on
   !> f(w) = w*exp(w) - x.  Returns zero for x <= 0, which lets callers detect
   !> a degenerate budget by testing the result.
   pure function lambert_w0(x) result(w)

      !> Argument, expected non-negative.
      real(dp), intent(in) :: x

      !> Solution of w*exp(w) = x on the principal branch.
      real(dp) :: w

      real(dp) :: exp_w, residual, w_plus_one, correction
      integer  :: iteration

      if (x <= 0.0_dp) then
         w = 0.0_dp
         return
      end if

      ! log(1 + x) lies above W0(x) for every x > 0, so the iteration
      ! approaches the root from one side and needs no bracketing stage.
      w = log(1.0_dp + x)

      do iteration = 1, max_iterations
         exp_w = exp(w)
         residual = w*exp_w - x
         w_plus_one = w + 1.0_dp
         correction = residual/(exp_w*w_plus_one - 0.5_dp*(w + 2.0_dp)*residual/w_plus_one)
         w = w - correction
         if (abs(correction) <= converged*(1.0_dp + abs(w))) exit
      end do

   end function lambert_w0

end module ewald_lambert_w

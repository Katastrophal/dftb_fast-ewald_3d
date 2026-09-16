module ewald_truncation
   !> Turning a requested accuracy into the two cutoffs that truncate an Ewald
   !> sum: the real-space radius r_cut and the Fourier-space radius k_cut.
   !>
   !> Two levels of estimate live here.  direct_cutoffs keeps only the Gaussian
   !> factor of each truncated tail and drops every prefactor; it is what the
   !> direct O(N^2) reference uses, where being conservative costs nothing
   !> that matters.  splitting_from_real_budget and cutoff_from_fourier_budget
   !> are root-mean-square estimates for a system of randomly placed charges
   !> and do keep the prefactors, which is what lets the fast methods place
   !> their cutoffs tightly.  Each is inverted for one unknown at a time and
   !> lands on a Lambert W value.
   !>
   !> Both refined estimates take the share of the requested accuracy that this
   !> truncation may spend as an argument rather than deriving it: the two fast
   !> methods do not split the request the same way.  See their parameter
   !> modules.
   !>
   !> They were also derived for a cubic, three-dimensionally periodic box.  A
   !> general cell has no edge length, so the caller passes the edge of the cube
   !> it considers equivalent to its own geometry.
   use ewald_constants, only: dp, pi
   use ewald_lambert_w, only: lambert_w0
   implicit none

   private
   public :: direct_cutoffs, splitting_from_real_budget, cutoff_from_fourier_budget

contains

   !> Elementary cutoffs for a prescribed splitting parameter: solve
   !> exp(-(alpha*r_cut)^2) = tol and exp(-(k_cut/(2*alpha))^2) = tol, which
   !> places each cutoff where the largest discarded term has the requested
   !> size.  Dropping the prefactors makes the result a target rather than a
   !> bound, but also makes it independent of the configuration.
   pure subroutine direct_cutoffs(alpha, tol, r_cut, k_cut)
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: tol
      real(dp), intent(out) :: r_cut
      real(dp), intent(out) :: k_cut

      r_cut = sqrt(-log(tol))/alpha
      k_cut = 2.0_dp*alpha*sqrt(-log(tol))

   end subroutine direct_cutoffs

   !> Splitting parameter that keeps the real-space truncation error inside a
   !> given budget at a prescribed cutoff radius.
   !>
   !> The estimate decays like exp(-(alpha*r_cut)^2) with a prefactor that
   !> grows with the charges and falls with the cell size; setting it equal to
   !> the budget and solving for (alpha*r_cut)^2 gives a Lambert W value.
   pure function splitting_from_real_budget(r_cut, chargeSquareSum, nParticle, &
                                            pseudoVolume, budget) result(alpha)
      real(dp), intent(in) :: r_cut
      real(dp), intent(in) :: chargeSquareSum
      integer, intent(in) :: nParticle
      real(dp), intent(in) :: pseudoVolume
      real(dp), intent(in) :: budget
      real(dp) :: alpha

      real(dp) :: argument, wValue

      argument = (2.0_dp/budget)*chargeSquareSum &
                 *sqrt(r_cut/(real(nParticle, dp)*pseudoVolume))

      if (argument > 0.0_dp) then
         wValue = lambert_w0(argument)
      else
         ! Degenerate budget or configuration: fall back on the elementary
         ! estimate, which asks only that the Gaussian factor sit at the budget.
         wValue = -log(budget)
      end if
      if (wValue <= 0.0_dp) wValue = -log(budget)

      alpha = sqrt(wValue)/r_cut

   end function splitting_from_real_budget

   !> Fourier-space cutoff that keeps the mode truncation error inside a given
   !> budget at a prescribed splitting parameter.  Same construction as the
   !> real-space inversion.  A caller working in mode indices converts the
   !> result with m = k_cut*|a_i|/(2*pi) along axis i.
   pure function cutoff_from_fourier_budget(chargeSquareSum, nParticle, boxEdge, &
                                            alpha, budget) result(k_cut)
      real(dp), intent(in) :: chargeSquareSum
      integer, intent(in) :: nParticle
      real(dp), intent(in) :: boxEdge
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: budget
      real(dp) :: k_cut

      real(dp) :: argument, wValue

      argument = (4.0_dp/(3.0_dp*boxEdge*boxEdge)) &
                 *(2.0_dp/(real(nParticle, dp)*alpha*pi))**(2.0_dp/3.0_dp) &
                 *(2.0_dp*chargeSquareSum/budget)**(4.0_dp/3.0_dp)

      wValue = lambert_w0(argument)
      if (wValue <= 0.0_dp) wValue = -log(budget)

      k_cut = sqrt(3.0_dp)*alpha*sqrt(wValue)

   end function cutoff_from_fourier_budget

end module ewald_truncation

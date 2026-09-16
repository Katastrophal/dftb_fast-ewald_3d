module ewald_alpha_search
   !> Choosing the Ewald splitting parameter by balancing the two truncated
   !> sums against each other.
   !>
   !> The total energy does not depend on the splitting parameter, so it cannot
   !> be chosen by minimising an error.  It is chosen instead by making the two
   !> branches converge at the same rate, taken concretely as: the drop of the
   !> summand between two successive shells is the same in both.  That
   !> condition is monotone in alpha, so it is solved by bracketing and
   !> bisection.  The direct 3d reference uses this search and differs only in
   !> how a single Fourier mode contributes.
   !>
   !> A tiny or very anisotropic cell can leave the condition without a sign
   !> change, in which case there is nothing to bisect.  The search then returns
   !> a non-zero status and the caller substitutes its own estimate; only the
   !> split of the work between the branches suffers, not the energy.
   use ewald_constants, only: dp, pi
   implicit none

   private
   public :: real_shell_term, balanced_alpha
   public :: alphaSearchStatus

   !> Outcome codes of the search.  Anything but `converged` means the caller
   !> has to supply its own splitting parameter.
   type :: TAlphaSearchStatus

      !> A bracket was found and bisected to within the tolerance.
      integer :: converged = 0

      !> The upward scan for the lower end of the bracket ran out of range.
      integer :: noLowerBracket = 1

      !> The condition never became negative, so the scan never started.
      integer :: alreadyBalanced = 2

      !> The upward scan for the upper end of the bracket ran out of range.
      integer :: noUpperBracket = 3

      !> Bisection did not reach the tolerance within the iteration cap.
      integer :: notConverged = 4

   end type TAlphaSearchStatus

   !> The status codes, so that callers read alphaSearchStatus%converged rather
   !> than a bare integer.
   type(TAlphaSearchStatus), parameter :: alphaSearchStatus = TAlphaSearchStatus()

   !> Start of the bracketing scan, deliberately far below any physical value
   !> so that the scan approaches the bracket from the unbalanced side.
   real(dp), parameter :: alphaStart = 1.0e-8_dp

   !> Cap on the number of bisection steps.
   integer, parameter :: maxBisection = 100

contains

   !> Magnitude of a single real-space pair term at distance r.  Common to both
   !> geometries, the real-space branch having the same form whatever the
   !> periodicity.
   pure function real_shell_term(r, alpha) result(termValue)

      !> Distance of the shell.
      real(dp), intent(in) :: r

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Magnitude of the term.
      real(dp) :: termValue

      termValue = erfc(alpha*r)/r

   end function real_shell_term

   !> Splitting parameter at which the two branches fall off at the same rate.
   !>
   !> The balance condition compares the drop of the Fourier summand between
   !> the fourth and fifth reciprocal shells with that of the real-space
   !> summand between the second and third real shells.  It is negative while
   !> the Fourier sum converges faster and positive once the real-space sum
   !> does, so a sign change brackets the balance point.
   subroutine balanced_alpha(minG, minR, cellMeasure, tolerance, &
                             alpha, status)

      !> Length of the shortest physical reciprocal lattice vector, including
      !> the factor of 2*pi.
      real(dp), intent(in) :: minG

      !> Length of the shortest real-space lattice vector.
      real(dp), intent(in) :: minR

      !> Volume of the 3d-periodic cell.
      real(dp), intent(in) :: cellMeasure

      !> How closely the condition has to be met.  Callers pass their accuracy
      !> request: there is no point balancing more finely than the accuracy the
      !> branches are truncated at.
      real(dp), intent(in) :: tolerance

      !> The balanced splitting parameter, meaningful only if status is
      !> alphaSearchStatus%converged.
      real(dp), intent(out) :: alpha

      !> Outcome of the search.
      integer, intent(out) :: status

      real(dp) :: alphaLeft, alphaRight, imbalance
      integer  :: iteration

      status = alphaSearchStatus%converged
      alpha = alphaStart
      imbalance = balance_condition(alpha)

      ! Scan upwards until the Fourier branch stops being the faster one.
      do while (imbalance < -tolerance .and. alpha <= huge(1.0_dp))
         alpha = 2.0_dp*alpha
         imbalance = balance_condition(alpha)
      end do

      if (alpha > huge(1.0_dp)) then
         status = alphaSearchStatus%noLowerBracket
      else if (alpha == alphaStart) then
         ! The condition was already non-negative at an unphysically small
         ! alpha, so there is no sign change to bracket.
         status = alphaSearchStatus%alreadyBalanced
      end if

      ! Keep scanning until the condition is clearly positive; the previous
      ! value then bounds the balance point from below.
      if (status == alphaSearchStatus%converged) then
         alphaLeft = 0.5_dp*alpha
         do while (imbalance < tolerance .and. alpha <= huge(1.0_dp))
            alpha = 2.0_dp*alpha
            imbalance = balance_condition(alpha)
         end do
         if (alpha > huge(1.0_dp)) status = alphaSearchStatus%noUpperBracket
      end if

      if (status == alphaSearchStatus%converged) then
         alphaRight = alpha
         alpha = (alphaLeft + alphaRight)/2.0_dp
         iteration = 0
         imbalance = balance_condition(alpha)
         do while (abs(imbalance) > tolerance .and. iteration <= maxBisection)
            if (imbalance < 0.0_dp) then
               alphaLeft = alpha
            else
               alphaRight = alpha
            end if
            alpha = (alphaLeft + alphaRight)/2.0_dp
            imbalance = balance_condition(alpha)
            iteration = iteration + 1
         end do
         if (iteration > maxBisection) status = alphaSearchStatus%notConverged
      end if

   contains

      !> How much faster the real-space sum converges than the Fourier sum at
      !> the given splitting: negative while the Fourier branch is ahead, zero
      !> at the balance point.
      function balance_condition(alphaTrial) result(imbalanceValue)
         real(dp), intent(in) :: alphaTrial
         real(dp) :: imbalanceValue

         real(dp) :: g4, g5   ! wavenumbers of the fourth and fifth shell
         real(dp) :: f4, f5   ! the weight the Fourier branch gives each

         g4 = 4.0_dp*minG
         g5 = 5.0_dp*minG

         f4 = (4.0_dp*pi/cellMeasure)*exp(-(g4**2)/(4.0_dp*alphaTrial**2))/g4**2
         f5 = (4.0_dp*pi/cellMeasure)*exp(-(g5**2)/(4.0_dp*alphaTrial**2))/g5**2

         imbalanceValue = (f4 - f5) &
                          - (real_shell_term(2.0_dp*minR, alphaTrial) &
                             - real_shell_term(3.0_dp*minR, alphaTrial))
      end function balance_condition

   end subroutine balanced_alpha

end module ewald_alpha_search

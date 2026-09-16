module ewald_slab_kernel
   !> The regularised kernel that lets a two-dimensionally periodic Ewald sum be
   !> evaluated by a three-dimensional FFT.
   !>
   !> In the fully periodic case the Fourier branch attaches a scalar weight to
   !> each mode and the sum factorises into squared structure factors, which is
   !> what an FFT computes.  In two dimensions the weight still depends on the
   !> separation z along the open direction,
   !>
   !>   K(kappa, z) = (pi/(kappa*A)) * Theta_+(kappa, z)
   !>
   !> for a non-zero in-plane wavenumber kappa, with a different expression at
   !> kappa = 0, so the sum no longer factorises and the direct reference is
   !> left evaluating it pair by pair.
   !>
   !> Factorisation is restored by transforming along z as well, which requires
   !> the kernel to be periodic in z, which it is not: it grows away from the
   !> slab.  It is therefore replaced by a function that agrees with it exactly
   !> wherever the physics needs it and is made periodic elsewhere.  A period h
   !> is chosen from the computed thickness so that the region where
   !> the two agree covers every separation that actually occurs; outside it, in
   !> a boundary layer of relative width eps at the edge of the period, the
   !> kernel is replaced by the polynomial that joins its two ends smoothly.
   !>
   !> The joint is a two-point Hermite interpolation matching the kernel's value
   !> and its first p-1 derivatives at both ends of the layer.  That smoothness
   !> is what makes the Fourier coefficients decay quickly and hence what keeps
   !> the mode set along z small enough to be worth transforming.  Two
   !> consequences reach the callers: the mode count along the open direction
   !> has to resolve the boundary layer and is therefore bounded from below by a
   !> multiple of p, and the regularised kernel is only p-1 times differentiable
   !> across the joint, so differentiating it for a force costs one order of
   !> smoothness the energy does not pay.
   !>
   !> Written out, each theta function contains a growing exponential multiplied
   !> by a decaying complementary error function, both of which overflow or
   !> underflow long before their product does.  They are therefore evaluated as
   !> erfc(x) = exp(-x^2)*erfcx(x), which pulls a single common Gaussian factor
   !> out in front.
   use ewald_constants, only: dp, pi, sqrt_pi
   use fft_backend, only: fft_1d_inplace
   implicit none

   private
   public :: regularised_kernel, kernel_fourier_coefficients
   public :: theta_plus, theta_minus, theta_zero

   !> In-plane wavenumber below which a mode is treated as the zero mode.
   real(dp), parameter :: minWavenumber = 1.0e-14_dp

contains

   ! =====================================================================
   !  The kernel and its building blocks
   ! =====================================================================

   !> The kernel at in-plane wavenumber kappa and separation r along the open
   !> direction, inside the region where no regularisation is applied.
   !>
   !> The zero mode is a separate expression rather than a limit: the
   !> pi/(kappa*A) prefactor diverges as kappa goes to zero while Theta_+
   !> vanishes, and the finite product is the uniform-sheet term carried by
   !> Theta_0.
   pure function regularised_kernel(kappa, r, area, alpha) result(value)

      !> In-plane wavenumber, |2 pi (m1 b1 + m2 b2)|.
      real(dp), intent(in) :: kappa

      !> Separation along the open direction.
      real(dp), intent(in) :: r

      !> Cell area.
      real(dp), intent(in) :: area

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Value of the kernel.
      real(dp) :: value

      if (kappa < minWavenumber) then
         value = -(2.0_dp*sqrt_pi/area)*theta_zero(alpha, r)
      else
         value = (pi/(kappa*area))*theta_plus(kappa, alpha, r)
      end if

   end function regularised_kernel

   !> Theta_+, the even combination
   !>
   !>   e^{kappa r} erfc(kappa/(2 alpha) + alpha r)
   !>     + e^{-kappa r} erfc(kappa/(2 alpha) - alpha r),
   !>
   !> evaluated in the scaled form described in the module header.
   pure function theta_plus(kappa, alpha, r) result(value)

      !> In-plane wavenumber.
      real(dp), intent(in) :: kappa

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Separation along the open direction.
      real(dp), intent(in) :: r

      !> Value of Theta_+.
      real(dp) :: value

      real(dp) :: gaussianFactor   ! the factor common to both terms

      gaussianFactor = exp(-(kappa/(2.0_dp*alpha))**2 - (alpha*r)**2)
      value = gaussianFactor*(erfc_scaled(kappa/(2.0_dp*alpha) + alpha*r) &
                              + erfc_scaled(kappa/(2.0_dp*alpha) - alpha*r))

   end function theta_plus

   !> Theta_-, the odd combination with a minus sign between the two terms.
   !>
   !> It is what the derivative of Theta_+ with respect to the separation
   !> reduces to: the Gaussian pieces the product rule generates cancel between
   !> the two terms exactly, leaving only the change of sign.
   pure function theta_minus(kappa, alpha, r) result(value)

      !> In-plane wavenumber.
      real(dp), intent(in) :: kappa

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Separation along the open direction.
      real(dp), intent(in) :: r

      !> Value of Theta_-.
      real(dp) :: value

      real(dp) :: gaussianFactor

      gaussianFactor = exp(-(kappa/(2.0_dp*alpha))**2 - (alpha*r)**2)
      value = gaussianFactor*(erfc_scaled(kappa/(2.0_dp*alpha) + alpha*r) &
                              - erfc_scaled(kappa/(2.0_dp*alpha) - alpha*r))

   end function theta_minus

   !> Theta_0, the zero-wavenumber limit,
   !>
   !>   e^{-alpha^2 r^2}/alpha + sqrt(pi) r erf(alpha r).
   !>
   !> It is the uniform-sheet term: what a charge feels from the others once
   !> their in-plane structure has been averaged away.
   pure function theta_zero(alpha, r) result(value)

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Separation along the open direction.
      real(dp), intent(in) :: r

      !> Value of Theta_0.
      real(dp) :: value

      value = exp(-(alpha*r)**2)/alpha + sqrt_pi*r*erf(alpha*r)

   end function theta_zero

   ! =====================================================================
   !  Fourier coefficients of the regularised, periodic kernel
   ! =====================================================================

   !> Discrete Fourier coefficients of the kernel at one in-plane wavenumber,
   !> made periodic over the period h.
   !>
   !> The kernel is sampled on a uniform grid over one period, taking its true
   !> value inside the physical region and the Hermite polynomial inside the
   !> boundary layer.  The coefficients come out real, the sampled function
   !> being even in the separation, and are returned in FFT bin order.
   !>
   !> Called once per in-plane mode, which is what makes the plan caching in the
   !> FFT backend's one-dimensional transform worth having.
   subroutine kernel_fourier_coefficients(kappa, area, alpha, period, eps, p, nModes, &
                                          coefficients)

      !> In-plane wavenumber of the mode.
      real(dp), intent(in) :: kappa

      !> Cell area.
      real(dp), intent(in) :: area

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Period imposed along the open direction.
      real(dp), intent(in) :: period

      !> Relative width of the boundary layer, as a fraction of the period.
      real(dp), intent(in) :: eps

      !> Smoothness of the regularisation: the polynomial matches the kernel's
      !> value and its first p-1 derivatives at both ends of the layer.
      integer, intent(in) :: p

      !> Number of samples, and hence of coefficients, along the period.
      integer, intent(in) :: nModes

      !> The coefficients, in FFT bin order.
      real(dp), intent(out) :: coefficients(0:nModes - 1)

      real(dp) :: physicalHalfWidth      ! separations up to here are untouched
      real(dp) :: layerCentre            ! centre of the boundary layer
      real(dp) :: layerHalfWidth         ! half width of the boundary layer
      real(dp) :: separation             ! sample point, reduced to one period
      real(dp) :: layerCoordinate        ! the sample point as the polynomial sees it
      real(dp) :: upperDerivatives(0:p - 1)   ! kernel and derivatives at +halfWidth
      real(dp) :: lowerDerivatives(0:p - 1)   ! kernel and derivatives at -halfWidth
      complex(dp), allocatable :: samples(:)
      integer :: iSample
      logical :: derivativesReady

      physicalHalfWidth = (0.5_dp - eps)*period
      layerCentre = 0.5_dp*period
      layerHalfWidth = eps*period

      derivativesReady = .false.
      allocate (samples(0:nModes - 1))

      do iSample = 0, nModes - 1
         separation = real(iSample, dp)*period/real(nModes, dp)
         ! Reduce to the symmetric interval, so that the physical region sits
         ! around zero and the boundary layer around the period's edge.
         separation = separation - period*nint(separation/period)

         if (abs(separation) <= physicalHalfWidth) then
            samples(iSample) = cmplx(regularised_kernel(kappa, separation, area, alpha), &
                                     0.0_dp, dp)
         else
            ! The end values are the same for every sample in the layer, so
            ! they are computed once and only if the layer is reached at all.
            if (.not. derivativesReady) then
               call boundary_derivatives(kappa, area, alpha, physicalHalfWidth, p, &
                                         upperDerivatives, lowerDerivatives)
               derivativesReady = .true.
            end if
            ! Samples on the negative side belong to the same layer, one period
            ! up; shifting them there is what makes the polynomial single valued.
            if (separation > 0.0_dp) then
               layerCoordinate = separation
            else
               layerCoordinate = separation + period
            end if
            samples(iSample) = cmplx(hermite_interpolant(layerCoordinate, layerCentre, &
                                                         layerHalfWidth, p, &
                                                         upperDerivatives, &
                                                         lowerDerivatives), 0.0_dp, dp)
         end if
      end do

      call fft_1d_inplace(samples, nModes, -1)
      coefficients = real(samples, dp)/real(nModes, dp)

   end subroutine kernel_fourier_coefficients

   !> Value and first p-1 derivatives of the kernel at the two ends of the
   !> physical region, which is what the interpolating polynomial has to match.
   subroutine boundary_derivatives(kappa, area, alpha, halfWidth, p, &
                                   upperDerivatives, lowerDerivatives)

      !> In-plane wavenumber of the mode.
      real(dp), intent(in) :: kappa

      !> Cell area.
      real(dp), intent(in) :: area

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Half width of the physical region: the ends are at plus and minus it.
      real(dp), intent(in) :: halfWidth

      !> Number of matched derivatives, counting the value itself.
      integer, intent(in) :: p

      !> Derivatives of order 0..p-1 at the upper end.
      real(dp), intent(out) :: upperDerivatives(0:p - 1)

      !> Derivatives of order 0..p-1 at the lower end.
      real(dp), intent(out) :: lowerDerivatives(0:p - 1)

      if (kappa < minWavenumber) then
         call theta_zero_derivatives(alpha, halfWidth, p, upperDerivatives)
         call theta_zero_derivatives(alpha, -halfWidth, p, lowerDerivatives)
         upperDerivatives = -(2.0_dp*sqrt_pi/area)*upperDerivatives
         lowerDerivatives = -(2.0_dp*sqrt_pi/area)*lowerDerivatives
      else
         call theta_plus_derivatives(kappa, alpha, halfWidth, p, upperDerivatives)
         call theta_plus_derivatives(kappa, alpha, -halfWidth, p, lowerDerivatives)
         upperDerivatives = (pi/(kappa*area))*upperDerivatives
         lowerDerivatives = (pi/(kappa*area))*lowerDerivatives
      end if

   end subroutine boundary_derivatives

   ! =====================================================================
   !  Derivatives of the theta functions
   ! =====================================================================

   !> Derivatives of Theta_+ with respect to the separation, orders 0..p-1.
   !>
   !> Differentiating twice returns Theta_+ multiplied by the squared wavenumber
   !> plus a Gaussian remainder, which gives a two-term recurrence: the whole
   !> ladder follows from the first two derivatives and the derivatives of a
   !> Gaussian.  Evaluating each order from its closed form would be slower and
   !> less stable.
   pure subroutine theta_plus_derivatives(kappa, alpha, r, p, derivatives)

      !> In-plane wavenumber.
      real(dp), intent(in) :: kappa

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Separation at which to differentiate.
      real(dp), intent(in) :: r

      !> Highest order needed, plus one.
      integer, intent(in) :: p

      !> Derivatives of orders 0..p-1.
      real(dp), intent(out) :: derivatives(0:p - 1)

      real(dp) :: remainderScale                 ! prefactor of the Gaussian remainder
      real(dp) :: gaussianDerivative(0:p - 1)    ! derivatives of exp(-alpha^2 r^2)
      integer  :: n

      call gaussian_derivatives(alpha, r, p, gaussianDerivative)
      remainderScale = (4.0_dp*alpha/sqrt_pi)*exp(-(kappa/(2.0_dp*alpha))**2)

      if (p >= 1) derivatives(0) = theta_plus(kappa, alpha, r)
      if (p >= 2) derivatives(1) = kappa*theta_minus(kappa, alpha, r)
      do n = 2, p - 1
         derivatives(n) = kappa**2*derivatives(n - 2) &
                          - kappa*remainderScale*gaussianDerivative(n - 2)
      end do

   end subroutine theta_plus_derivatives

   !> Derivatives of Theta_0 with respect to the separation, orders 0..p-1.  The
   !> recurrence of theta_plus_derivatives with the wavenumber set to zero: its
   !> first term drops out and every order above the first is a Gaussian
   !> derivative.
   pure subroutine theta_zero_derivatives(alpha, r, p, derivatives)

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Separation at which to differentiate.
      real(dp), intent(in) :: r

      !> Highest order needed, plus one.
      integer, intent(in) :: p

      !> Derivatives of orders 0..p-1.
      real(dp), intent(out) :: derivatives(0:p - 1)

      real(dp) :: gaussianDerivative(0:p - 1)
      integer  :: n

      call gaussian_derivatives(alpha, r, p, gaussianDerivative)

      if (p >= 1) derivatives(0) = theta_zero(alpha, r)
      if (p >= 2) derivatives(1) = sqrt_pi*erf(alpha*r)
      do n = 2, p - 1
         derivatives(n) = 2.0_dp*alpha*gaussianDerivative(n - 2)
      end do

   end subroutine theta_zero_derivatives

   !> Derivatives of exp(-alpha^2 r^2), orders 0..p-1.  Each is the Gaussian
   !> multiplied by a Hermite polynomial of the scaled argument, up to a sign
   !> and a power of alpha.  The polynomials are built by their own three-term
   !> recurrence rather than evaluated from their coefficients, which would lose
   !> precision at the orders used here.
   pure subroutine gaussian_derivatives(alpha, r, p, derivatives)

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Point at which to differentiate.
      real(dp), intent(in) :: r

      !> Highest order needed, plus one.
      integer, intent(in) :: p

      !> Derivatives of orders 0..p-1.
      real(dp), intent(out) :: derivatives(0:p - 1)

      real(dp) :: scaledArgument                  ! alpha * r
      real(dp) :: hermite(0:max(p - 1, 1))        ! Hermite polynomials at that point
      real(dp) :: gaussian                        ! exp(-(alpha r)^2)
      integer  :: n

      scaledArgument = alpha*r
      gaussian = exp(-scaledArgument*scaledArgument)

      hermite(0) = 1.0_dp
      if (p - 1 >= 1) hermite(1) = 2.0_dp*scaledArgument
      do n = 1, p - 2
         hermite(n + 1) = 2.0_dp*scaledArgument*hermite(n) - 2.0_dp*real(n, dp)*hermite(n - 1)
      end do

      do n = 0, p - 1
         derivatives(n) = (alpha**n)*((-1.0_dp)**n)*hermite(n)*gaussian
      end do

   end subroutine gaussian_derivatives

   ! =====================================================================
   !  Two-point Hermite interpolation
   ! =====================================================================

   !> The interpolating polynomial evaluated at one point of the boundary layer.
   !>
   !> It is written in a basis whose members vanish, together with all their
   !> relevant derivatives, at one end of the layer while reproducing one
   !> prescribed derivative at the other, so the polynomial is simply the
   !> prescribed derivatives weighted by that basis, with no linear system to
   !> solve.  The powers of the half width convert between the derivatives of
   !> the kernel and those of the scaled variable the basis is written in.
   pure function hermite_interpolant(x, centre, halfWidth, p, &
                                     upperDerivatives, lowerDerivatives) result(value)

      !> Point in the boundary layer at which to evaluate.
      real(dp), intent(in) :: x

      !> Centre of the boundary layer.
      real(dp), intent(in) :: centre

      !> Half width of the boundary layer.
      real(dp), intent(in) :: halfWidth

      !> Number of matched derivatives at each end, counting the value.
      integer, intent(in) :: p

      !> Derivatives of orders 0..p-1 prescribed at the upper end.
      real(dp), intent(in) :: upperDerivatives(0:p - 1)

      !> Derivatives of orders 0..p-1 prescribed at the lower end.
      real(dp), intent(in) :: lowerDerivatives(0:p - 1)

      !> Value of the interpolating polynomial.
      real(dp) :: value

      real(dp) :: scaled     ! the point mapped onto [-1, 1]
      integer  :: j          ! order of the derivative being matched

      scaled = (x - centre)/halfWidth
      value = 0.0_dp
      do j = 0, p - 1
         value = value + hermite_basis(p, j, scaled)*(halfWidth**j)*upperDerivatives(j) &
                 + hermite_basis(p, j, -scaled)*((-halfWidth)**j)*lowerDerivatives(j)
      end do

   end function hermite_interpolant

   !> One member of the interpolation basis, of degree 2p-1 in the scaled
   !> variable.  It has a zero of order p at one end of the interval, which is
   !> what makes it invisible to every condition imposed there.
   pure function hermite_basis(p, j, y) result(value)

      !> Number of matched derivatives at each end.
      integer, intent(in) :: p

      !> Order of the derivative this member reproduces.
      integer, intent(in) :: j

      !> Point in the scaled interval [-1, 1].
      real(dp), intent(in) :: y

      !> Value of the basis polynomial.
      real(dp) :: value

      real(dp) :: term
      integer  :: k

      value = 0.0_dp
      do k = 0, p - 1 - j
         term = binomial(p - 1 + k, k)/(factorial(j)*(2.0_dp**p)*(2.0_dp**k))
         term = term*(1.0_dp - y)**p*(1.0_dp + y)**(k + j)
         value = value + term
      end do

   end function hermite_basis

   !> Binomial coefficient, built multiplicatively so that no intermediate
   !> factorial has to be represented.
   pure function binomial(n, k) result(value)

      !> Upper index.
      integer, intent(in) :: n

      !> Lower index.
      integer, intent(in) :: k

      !> The coefficient, as a real number.
      real(dp) :: value

      integer :: i

      value = 1.0_dp
      do i = 0, k - 1
         value = value*real(n - i, dp)/real(i + 1, dp)
      end do

   end function binomial

   !> Factorial, as a real number.
   pure function factorial(n) result(value)

      !> Argument, expected non-negative.
      integer, intent(in) :: n

      !> The factorial.
      real(dp) :: value

      integer :: i

      value = 1.0_dp
      do i = 2, n
         value = value*real(i, dp)
      end do

   end function factorial

end module ewald_slab_kernel

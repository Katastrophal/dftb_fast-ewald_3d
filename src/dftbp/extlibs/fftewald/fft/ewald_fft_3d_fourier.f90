module ewald_fft_3d_fourier
   !> Fourier branch of the fast three-dimensional potential/force method.
   !>
   !> A non-uniform transform obtains all signed structure factors in
   !> O(N + M log M). The weighted coefficients are interpolated to the
   !> particles for the potential and differentiated for the force.
   use ewald_constants, only: dp, pi
   use ewald_validation, only: minWavenumberSquared
   use nfft, only: adjoint_nfft_3d, forward_nfft_3d, &
                   mode_of_bin, oversampling, window_shape
   implicit none

   private
   public :: long_range_potential_force

contains

   !> Fourier contribution to the per-atom potential and force.
   !>
   !> One adjoint transform obtains the structure factors, one forward
   !> transform of the weighted coefficients gives the potential, and one
   !> forward transform per force component gives the force; the force
   !> transforms are skipped when the caller asks for no force.
   !> Differentiating the potential turns each mode's phase into a factor of
   !> the wavevector, which is why the force comes from transforming the same
   !> coefficients multiplied by a component of k.  The results are real up to
   !> round-off, the conjugate mode pairs cancelling the imaginary parts.
   subroutine long_range_potential_force(positions, charges, nParticle, recVecs, volume, &
                                         alpha, nModes, windowCutoff, pot, force)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(3, 3)

      !> Cell volume.
      real(dp), intent(in) :: volume

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Number of retained modes along each axis.
      integer, intent(in) :: nModes(3)

      !> Half-width of the transform's window stencil.
      integer, intent(in) :: windowCutoff

      !> Fourier part of the potential at each charge.
      real(dp), intent(out) :: pot(:)

      !> Fourier part of the force on each charge, shape (3, nParticle).
      !> Omitting it skips the three extra transforms the force costs.
      real(dp), intent(out), optional :: force(:, :)

      real(dp), allocatable :: t1(:), t2(:), t3(:)
      complex(dp), allocatable :: coefficients(:, :, :)   ! Shat, then b(k)*Shat
      complex(dp), allocatable :: derivative(:, :, :)     ! one force component's coefficients
      complex(dp), allocatable :: atSites(:)              ! a transform evaluated at the charges
      integer  :: i
      integer  :: a1, a2, a3          ! mode bins
      integer  :: k1, k2, k3          ! signed modes
      integer  :: iComponent          ! force component being transformed
      real(dp) :: waveVector(3)
      real(dp) :: wavenumberSquared
      real(dp) :: modeWeight
      real(dp) :: prefactor
      real(dp) :: inverseFourAlphaSq

      allocate (t1(nParticle), t2(nParticle), t3(nParticle))
      do i = 1, nParticle
         t1(i) = dot_product(recVecs(1, :), positions(i, :))
         t2(i) = dot_product(recVecs(2, :), positions(i, :))
         t3(i) = dot_product(recVecs(3, :), positions(i, :))
      end do

      ! Signed structure factors are required for the potential and force, so
      ! the full complex transform is used.
      allocate (coefficients(0:nModes(1) - 1, 0:nModes(2) - 1, 0:nModes(3) - 1))
      call adjoint_nfft_3d(nParticle, t1, t2, t3, charges(1:nParticle), nModes(1), nModes(2), nModes(3), &
                           windowCutoff, coefficients)

      ! Weight the structure factors in place.  Overwriting rather than
      ! allocating a second array of this size matters: it is the largest
      ! complex array of the evaluation.
      prefactor = 4.0_dp*pi/volume
      inverseFourAlphaSq = 1.0_dp/(4.0_dp*alpha**2)
      do a3 = 0, nModes(3) - 1
         k3 = mode_of_bin(a3, nModes(3))
         do a2 = 0, nModes(2) - 1
            k2 = mode_of_bin(a2, nModes(2))
            do a1 = 0, nModes(1) - 1
               k1 = mode_of_bin(a1, nModes(1))
               if (k1 == 0 .and. k2 == 0 .and. k3 == 0) then
                  coefficients(a1, a2, a3) = (0.0_dp, 0.0_dp)
                  cycle
               end if

               waveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                       + real(k2, dp)*recVecs(2, :) &
                                       + real(k3, dp)*recVecs(3, :))
               wavenumberSquared = sum(waveVector**2)
               if (wavenumberSquared < minWavenumberSquared) then
                  coefficients(a1, a2, a3) = (0.0_dp, 0.0_dp)
                  cycle
               end if

               modeWeight = prefactor*exp(-wavenumberSquared*inverseFourAlphaSq) &
                            /wavenumberSquared
               coefficients(a1, a2, a3) = modeWeight*coefficients(a1, a2, a3)
            end do
         end do
      end do
      ! Potential: one forward transform of the weighted coefficients.
      allocate (atSites(nParticle))
      call forward_nfft_3d(nParticle, t1, t2, t3, nModes(1), nModes(2), nModes(3), &
                           windowCutoff, coefficients, atSites)
      pot(1:nParticle) = real(atSites(1:nParticle), dp)

      if (.not. present(force)) return

      ! Force: one forward transform per component, of the same coefficients
      ! multiplied by that component of the wavevector.
      allocate (derivative(0:nModes(1) - 1, 0:nModes(2) - 1, 0:nModes(3) - 1))
      do iComponent = 1, 3
         do a3 = 0, nModes(3) - 1
            k3 = mode_of_bin(a3, nModes(3))
            do a2 = 0, nModes(2) - 1
               k2 = mode_of_bin(a2, nModes(2))
               do a1 = 0, nModes(1) - 1
                  k1 = mode_of_bin(a1, nModes(1))
                  waveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                          + real(k2, dp)*recVecs(2, :) &
                                          + real(k3, dp)*recVecs(3, :))
                  derivative(a1, a2, a3) = cmplx(0.0_dp, waveVector(iComponent), dp) &
                                           *coefficients(a1, a2, a3)
               end do
            end do
         end do
         call forward_nfft_3d(nParticle, t1, t2, t3, nModes(1), nModes(2), nModes(3), &
                              windowCutoff, derivative, atSites)
         do i = 1, nParticle
            force(iComponent, i) = charges(i)*real(atSites(i), dp)
         end do
      end do

   end subroutine long_range_potential_force

end module ewald_fft_3d_fourier

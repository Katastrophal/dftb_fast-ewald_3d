module ewald_fft_3d_fourier
   !> Fourier branch of the fast three-dimensional method.
   !>
   !> A non-uniform transform obtains all structure factors in
   !> O(N + M log M). The Fourier contribution is
   !>
   !>   U = (1/2) sum_{k /= 0} b(k) |Shat(k)|^2 ,
   !>   b(k) = (4 pi / V) exp(-k^2/(4 alpha^2)) / k^2 ,
   !>
   !> which is the expression the direct reference evaluates; the two differ
   !> only in how Shat is obtained.
   !>
   !> The energy needs only the magnitudes of the structure factors and the
   !> spread grid is purely real, so its path spreads onto a padded real grid,
   !> runs a real-input FFT that stores only the non-redundant half of the
   !> spectrum, and folds the window deconvolution into the mode sum.  The
   !> potential and force need the signed structure factors, so they run the
   !> general complex transforms: one for the potential and one per force
   !> component.
   use ewald_constants, only: dp, pi
   use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
   use ewald_validation, only: minWavenumberSquared
   use fft_backend, only: fft_3d_real_to_complex
   use nfft, only: spread_charges_real_3d, adjoint_nfft_3d, forward_nfft_3d, &
                   mode_of_bin, oversampling, window_shape
   implicit none

   private
   public :: long_range_energy, long_range_potential_force

contains

   !> Fourier contribution to the energy.
   function long_range_energy(positions, charges, nParticle, recVecs, volume, alpha, &
                              nModes, windowCutoff, phaseTimes) result(energy)

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

      !> Wall time in seconds of the four phases, as
      !> [spread, transform, extract, mode sum].  The extract phase is folded
      !> into the mode sum here and is reported as zero.  For profiling only.
      real(dp), intent(out), optional :: phaseTimes(4)

      !> Fourier contribution to the total energy.
      real(dp) :: energy

      real(dp), allocatable :: t1(:), t2(:), t3(:)   ! fractional coordinates

      ! The one large allocation of the method: a padded real grid the charges
      ! are spread onto, overwritten in place by the real-input transform with
      ! the half spectrum aliased below.
      real(dp), allocatable, target :: realGrid(:, :, :)
      complex(dp), pointer :: spectrumBase(:, :, :)   ! complex view, default bounds
      complex(dp), pointer :: spectrum(:, :, :)       ! the same, indexed from zero

      real(dp), allocatable :: windowFactor1(:), windowFactor2(:), windowFactor3(:)
      integer  :: i
      integer  :: n1, n2, n3          ! fine grid extents
      integer  :: halfExtent1         ! stored extent of the half spectrum
      integer  :: a1, a2, a3          ! mode bins
      integer  :: k1, k2, k3          ! signed modes
      integer  :: j1, j2, j3          ! indices into the stored half spectrum
      real(dp) :: b                   ! window shape parameter
      real(dp) :: waveVector(3)
      real(dp) :: wavenumberSquared
      real(dp) :: modeWeight          ! b(k)
      real(dp) :: prefactor, inverseFourAlphaSq
      real(dp) :: deconvolution       ! product of the three window factors
      real(dp) :: magnitudeSquared    ! |Shat(k)|^2 before deconvolution
      integer(8) :: tick0, tick1, tick2, tick3, tickRate

      allocate (t1(nParticle), t2(nParticle), t3(nParticle))
      do i = 1, nParticle
         t1(i) = dot_product(recVecs(1, :), positions(i, :))
         t2(i) = dot_product(recVecs(2, :), positions(i, :))
         t3(i) = dot_product(recVecs(3, :), positions(i, :))
      end do

      n1 = oversampling*nModes(1)
      n2 = oversampling*nModes(2)
      n3 = oversampling*nModes(3)
      halfExtent1 = n1/2 + 1
      b = window_shape(windowCutoff)

      ! Window deconvolution factors, one table per axis.  They are real and
      ! even in the mode index, which is what lets the first axis be tabulated
      ! over the stored half alone: the table is indexed by |k1|.  The other
      ! two are tabulated by bin, so both signs are already covered.
      allocate (windowFactor1(0:nModes(1)/2), windowFactor2(0:nModes(2) - 1), &
                windowFactor3(0:nModes(3) - 1))
      do j1 = 0, nModes(1)/2
         windowFactor1(j1) = exp(b*(pi*real(j1, dp)/real(n1, dp))**2)
      end do
      do a2 = 0, nModes(2) - 1
         windowFactor2(a2) = exp(b*(pi*real(mode_of_bin(a2, nModes(2)), dp) &
                                    /real(n2, dp))**2)
      end do
      do a3 = 0, nModes(3) - 1
         windowFactor3(a3) = exp(b*(pi*real(mode_of_bin(a3, nModes(3)), dp) &
                                    /real(n3, dp))**2)
      end do

      ! The real grid and a complex view of the same storage.  Handing both to
      ! the transform is what makes it in place; the padding of the first
      ! dimension is the layout the real-input transform requires.
      allocate (realGrid(0:2*halfExtent1 - 1, 0:n2 - 1, 0:n3 - 1))
      realGrid = 0.0_dp
      call c_f_pointer(c_loc(realGrid), spectrumBase, [halfExtent1, n2, n3])
      spectrum(0:, 0:, 0:) => spectrumBase

      call system_clock(tick0, tickRate)
      call spread_charges_real_3d(nParticle, t1, t2, t3, charges(1:nParticle), n1, n2, n3, &
                                  windowCutoff, realGrid)
      call system_clock(tick1)
      call fft_3d_real_to_complex(realGrid, spectrum, n1, n2, n3)
      call system_clock(tick2)

      ! Mode sum, with the deconvolution folded in.  The loop visits the full
      ! mode set but reads the stored half: for a real grid the spectrum is
      ! Hermitian, and the weight and the deconvolution are both even in k, so
      ! mapping a mode onto its representative with non-negative k1 is exact.
      prefactor = 4.0_dp*pi/volume
      inverseFourAlphaSq = 1.0_dp/(4.0_dp*alpha**2)
      energy = 0.0_dp

      do a3 = 0, nModes(3) - 1
         k3 = mode_of_bin(a3, nModes(3))
         do a2 = 0, nModes(2) - 1
            k2 = mode_of_bin(a2, nModes(2))
            do a1 = 0, nModes(1) - 1
               k1 = mode_of_bin(a1, nModes(1))
               if (k1 == 0 .and. k2 == 0 .and. k3 == 0) cycle

               waveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                       + real(k2, dp)*recVecs(2, :) &
                                       + real(k3, dp)*recVecs(3, :))
               wavenumberSquared = sum(waveVector**2)
               if (wavenumberSquared < minWavenumberSquared) cycle

               modeWeight = prefactor*exp(-wavenumberSquared*inverseFourAlphaSq) &
                            /wavenumberSquared

               if (k1 >= 0) then
                  j1 = k1; j2 = modulo(k2, n2); j3 = modulo(k3, n3)
               else
                  j1 = -k1; j2 = modulo(-k2, n2); j3 = modulo(-k3, n3)
               end if

               deconvolution = windowFactor1(j1)*windowFactor2(a2)*windowFactor3(a3)
               magnitudeSquared = real(spectrum(j1, j2, j3))**2 &
                                  + aimag(spectrum(j1, j2, j3))**2
               energy = energy + modeWeight*deconvolution*deconvolution*magnitudeSquared
            end do
         end do
      end do
      energy = 0.5_dp*energy
      call system_clock(tick3)

      nullify (spectrum, spectrumBase)

      if (present(phaseTimes)) &
         phaseTimes = [real(tick1 - tick0, dp), real(tick2 - tick1, dp), &
                       0.0_dp, real(tick3 - tick2, dp)]/real(tickRate, dp)

   end function long_range_energy

   !> Fourier contribution to the per-atom potential and force, and to the
   !> energy as a by-product.
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
                                         alpha, nModes, windowCutoff, energy, pot, force)

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

      !> Fourier contribution to the total energy, identical to what
      !> long_range_energy returns for the same arguments.
      real(dp), intent(out) :: energy

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

      ! Signed structure factors.  The real-input shortcut of the energy path
      ! is unavailable here, so the full complex transform is used.
      allocate (coefficients(0:nModes(1) - 1, 0:nModes(2) - 1, 0:nModes(3) - 1))
      call adjoint_nfft_3d(nParticle, t1, t2, t3, charges(1:nParticle), nModes(1), nModes(2), nModes(3), &
                           windowCutoff, coefficients)

      ! Weight the structure factors in place, taking the energy from the same
      ! pass.  Overwriting rather than allocating a second array of this size
      ! matters: it is the largest complex array of the evaluation.
      prefactor = 4.0_dp*pi/volume
      inverseFourAlphaSq = 1.0_dp/(4.0_dp*alpha**2)
      energy = 0.0_dp

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
               energy = energy + modeWeight*(real(coefficients(a1, a2, a3))**2 &
                                             + aimag(coefficients(a1, a2, a3))**2)
               coefficients(a1, a2, a3) = modeWeight*coefficients(a1, a2, a3)
            end do
         end do
      end do
      energy = 0.5_dp*energy

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

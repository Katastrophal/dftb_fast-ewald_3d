module ewald_fft_2d_fourier
   !> Fourier branch of the fast two-dimensional method.
   !>
   !> Here the mode weight still depends on the separation along the open
   !> direction, so the sum factorises only after the kernel has been made
   !> periodic along that direction and transformed as well.  The regularisation
   !> that does so lives in ewald_slab_kernel; this module uses it.
   !>
   !> Two paths, one routine each for the energy and for the potential and
   !> force.  On a monolayer every charge lies in one plane, the only separation
   !> along the open direction is zero, the kernel is a scalar per in-plane mode
   !> and no regularisation is needed; the branch reduces to a two-dimensional
   !> transform with a single mode along the open direction.  On a slab a period
   !> is imposed along the open direction so that every separation that actually
   !> occurs falls in the region where the regularised kernel equals the true
   !> one, and the transform runs over three directions; the kernel's
   !> coefficients then depend on the mode along the open direction, so one
   !> small transform of the kernel is needed per in-plane mode.
   !>
   !> Which path applies, and the period and centre the regularisation was sized
   !> for, are settled by the parameter chain and arrive in TEwaldParameters2d.
   !> Nothing here re-measures the geometry, so the transform cannot disagree
   !> with the mode counts it was handed.
   !>
   !> The zero in-plane mode needs no special handling in either path: the
   !> kernel carries the uniform-sheet term there, so what is a separate
   !> quadratic double sum in the direct reference comes out of the same
   !> transform as every other mode.
   !>
   !> As in three dimensions, the energy needs only squared magnitudes and
   !> spreads onto a purely real grid, so it uses a real-input transform that
   !> stores half the spectrum and folds the window deconvolution into its own
   !> mode sum.  The potential and force paths need the signed structure factors
   !> and run the general complex transforms.
   use ewald_constants, only: dp, pi
   use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
   use fft_backend, only: fft_2d_real_to_complex, fft_3d_real_to_complex
   use nfft, only: spread_charges_real_2d, spread_charges_real_3d, &
                   adjoint_nfft_2d, adjoint_nfft_3d, &
                   forward_nfft_2d, forward_nfft_3d, &
                   mode_of_bin, oversampling, window_shape
   use ewald_slab_kernel, only: regularised_kernel, kernel_fourier_coefficients
   use ewald_fft_2d_parameters, only: TEwaldParameters2d
   implicit none

   private
   public :: long_range_energy
   public :: long_range_potential_force_monolayer, long_range_potential_force_slab

contains

   !> Fourier contribution to the energy, by whichever path the geometry
   !> selects.  The two paths share nothing but the in-plane window tables, so
   !> each has a routine of its own, exactly as the potential and force do.
   function long_range_energy(positions, charges, nParticle, recVecs, area, params, &
                              phaseTimes) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> In-plane crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(2, 3)

      !> Cell area.
      real(dp), intent(in) :: area

      !> Parameters of this evaluation, including the selected path and the
      !> period and centre for which the regularisation was sized.
      type(TEwaldParameters2d), intent(in) :: params

      !> Wall time in seconds of the four phases, as
      !> [spread, transform, extract, mode sum].  The extract phase is folded
      !> into the mode sum and reported as zero.  For profiling only.
      real(dp), intent(out), optional :: phaseTimes(4)

      !> Fourier contribution to the total energy.
      real(dp) :: energy

      if (params%monolayer) then
         energy = monolayer_energy(positions, charges, nParticle, recVecs, area, &
                                   params, phaseTimes)
      else
         energy = slab_energy(positions, charges, nParticle, recVecs, area, &
                              params, phaseTimes)
      end if

   end function long_range_energy

   !> Energy of a coplanar layer: one in-plane transform, and the kernel taken
   !> at zero separation along the open direction, which is exact here and is
   !> the reason no regularisation is needed.
   function monolayer_energy(positions, charges, nParticle, recVecs, area, params, &
                             phaseTimes) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> In-plane crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(2, 3)

      !> Cell area.
      real(dp), intent(in) :: area

      !> Parameters of this evaluation.
      type(TEwaldParameters2d), intent(in) :: params

      !> Wall time in seconds of the four phases.  For profiling only.
      real(dp), intent(out), optional :: phaseTimes(4)

      !> Fourier contribution to the total energy.
      real(dp) :: energy

      real(dp), allocatable :: t1(:), t2(:)   ! coordinates on the torus

      ! The one large allocation: a padded real grid the charges are spread
      ! onto, overwritten in place by the real-input transform with the half
      ! spectrum aliased below.
      real(dp), allocatable, target :: planeGrid(:, :)
      complex(dp), pointer :: spectrumBase(:, :), spectrum(:, :)

      real(dp), allocatable :: windowFactor1(:), windowFactor2(:)
      integer  :: i
      integer  :: n1, n2               ! fine grid extents
      integer  :: halfExtent1          ! stored extent of the half spectrum
      integer  :: a1, a2               ! mode bins
      integer  :: k1, k2               ! signed modes
      integer  :: j1, j2               ! indices into the stored half spectrum
      real(dp) :: b                    ! window shape parameter
      real(dp) :: waveVector(3)
      real(dp) :: wavenumber
      real(dp) :: modeWeight
      real(dp) :: deconvolution
      real(dp) :: magnitudeSquared
      integer(8) :: tick0, tick1, tick2, tick3, tickRate

      allocate (t1(nParticle), t2(nParticle))
      do i = 1, nParticle
         t1(i) = dot_product(recVecs(1, :), positions(i, :))
         t2(i) = dot_product(recVecs(2, :), positions(i, :))
      end do

      n1 = oversampling*params%nModes(1)
      n2 = oversampling*params%nModes(2)
      halfExtent1 = n1/2 + 1
      b = window_shape(params%windowCutoff)

      allocate (windowFactor1(0:params%nModes(1)/2), windowFactor2(0:params%nModes(2) - 1))
      call window_factors_half(params%nModes(1), n1, b, windowFactor1)
      call window_factors(params%nModes(2), n2, b, windowFactor2)

      allocate (planeGrid(0:2*halfExtent1 - 1, 0:n2 - 1))
      planeGrid = 0.0_dp
      call c_f_pointer(c_loc(planeGrid), spectrumBase, [halfExtent1, n2])
      spectrum(0:, 0:) => spectrumBase

      call system_clock(tick0, tickRate)
      call spread_charges_real_2d(nParticle, t1, t2, charges(1:nParticle), n1, n2, &
                                  params%windowCutoff, planeGrid)
      call system_clock(tick1)
      call fft_2d_real_to_complex(planeGrid, spectrum, n1, n2)
      call system_clock(tick2)

      ! The zero in-plane mode is kept: the kernel carries the uniform-sheet
      ! term there.  A mode with k1 < 0 has the same magnitude as its mirror
      ! image, the spectrum of a real grid being Hermitian, so reading the
      ! stored half at |k1| is exact.
      energy = 0.0_dp
      do a2 = 0, params%nModes(2) - 1
         k2 = mode_of_bin(a2, params%nModes(2))
         do a1 = 0, params%nModes(1) - 1
            k1 = mode_of_bin(a1, params%nModes(1))
            waveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                    + real(k2, dp)*recVecs(2, :))
            wavenumber = sqrt(sum(waveVector**2))
            modeWeight = regularised_kernel(wavenumber, 0.0_dp, area, params%alpha)

            if (k1 >= 0) then
               j1 = k1; j2 = modulo(k2, n2)
            else
               j1 = -k1; j2 = modulo(-k2, n2)
            end if

            deconvolution = windowFactor1(j1)*windowFactor2(a2)
            magnitudeSquared = real(spectrum(j1, j2))**2 + aimag(spectrum(j1, j2))**2
            energy = energy + modeWeight*deconvolution*deconvolution*magnitudeSquared
         end do
      end do
      energy = 0.5_dp*energy
      call system_clock(tick3)

      nullify (spectrum, spectrumBase)

      if (present(phaseTimes)) &
         phaseTimes = [real(tick1 - tick0, dp), real(tick2 - tick1, dp), &
                       0.0_dp, real(tick3 - tick2, dp)]/real(tickRate, dp)

   end function monolayer_energy

   !> Energy of a slab: the full three-dimensional transform, with the kernel
   !> made periodic over the period the parameter chain chose.  That period
   !> reaches exactly the computed thickness, which is the largest separation
   !> any pair can have, so nothing physical is affected by the regularisation.
   function slab_energy(positions, charges, nParticle, recVecs, area, params, &
                        phaseTimes) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> In-plane crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(2, 3)

      !> Cell area.
      real(dp), intent(in) :: area

      !> Parameters of this evaluation.
      type(TEwaldParameters2d), intent(in) :: params

      !> Wall time in seconds of the four phases.  For profiling only.
      real(dp), intent(out), optional :: phaseTimes(4)

      !> Fourier contribution to the total energy.
      real(dp) :: energy

      real(dp), allocatable :: t1(:), t2(:), t3(:)   ! coordinates on the torus

      ! The one large allocation, as on the monolayer path.
      real(dp), allocatable, target :: slabGrid(:, :, :)
      complex(dp), pointer :: spectrumBase(:, :, :), spectrum(:, :, :)

      real(dp), allocatable :: windowFactor1(:), windowFactor2(:), windowFactor3(:)
      real(dp), allocatable :: kernelCoefficients(:)
      integer  :: i
      integer  :: n1, n2, n3           ! fine grid extents
      integer  :: halfExtent1          ! stored extent of the half spectrum
      integer  :: a1, a2, a3           ! mode bins
      integer  :: k1, k2, k3           ! signed modes
      integer  :: j1, j2, j3           ! indices into the stored half spectrum
      integer  :: conjugateSign        ! folds the Hermitian mirror into the index map
      integer  :: normalModes          ! modes along the open direction, after the guard
      real(dp) :: b                    ! window shape parameter
      real(dp) :: waveVector(3)
      real(dp) :: wavenumber
      real(dp) :: deconvolution
      real(dp) :: magnitudeSquared
      integer(8) :: tick0, tick1, tick2, tick3, tickRate

      n1 = oversampling*params%nModes(1)
      n2 = oversampling*params%nModes(2)
      halfExtent1 = n1/2 + 1
      normalModes = max(params%nModes(3), 1)
      n3 = oversampling*normalModes
      b = window_shape(params%windowCutoff)

      allocate (t1(nParticle), t2(nParticle), t3(nParticle))
      do i = 1, nParticle
         t1(i) = dot_product(recVecs(1, :), positions(i, :))
         t2(i) = dot_product(recVecs(2, :), positions(i, :))
         t3(i) = (positions(i, 3) - params%centreZ)/params%period
      end do

      allocate (windowFactor1(0:params%nModes(1)/2), windowFactor2(0:params%nModes(2) - 1), &
                windowFactor3(0:normalModes - 1))
      call window_factors_half(params%nModes(1), n1, b, windowFactor1)
      call window_factors(params%nModes(2), n2, b, windowFactor2)
      call window_factors(normalModes, n3, b, windowFactor3)

      allocate (slabGrid(0:2*halfExtent1 - 1, 0:n2 - 1, 0:n3 - 1))
      slabGrid = 0.0_dp
      call c_f_pointer(c_loc(slabGrid), spectrumBase, [halfExtent1, n2, n3])
      spectrum(0:, 0:, 0:) => spectrumBase

      call system_clock(tick0, tickRate)
      call spread_charges_real_3d(nParticle, t1, t2, t3, charges(1:nParticle), n1, n2, n3, &
                                  params%windowCutoff, slabGrid)
      call system_clock(tick1)
      call fft_3d_real_to_complex(slabGrid, spectrum, n1, n2, n3)
      call system_clock(tick2)

      energy = 0.0_dp
      allocate (kernelCoefficients(0:normalModes - 1))
      do a2 = 0, params%nModes(2) - 1
         k2 = mode_of_bin(a2, params%nModes(2))
         do a1 = 0, params%nModes(1) - 1
            k1 = mode_of_bin(a1, params%nModes(1))
            waveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                    + real(k2, dp)*recVecs(2, :))
            wavenumber = sqrt(sum(waveVector**2))

            ! The kernel depends on the mode along the open direction, so its
            ! coefficients have to be transformed once per in-plane mode.
            call kernel_fourier_coefficients(wavenumber, area, params%alpha, &
                                             params%period, params%layerWidth, &
                                             params%smoothness, normalModes, &
                                             kernelCoefficients)

            ! Hermitian read: for k1 < 0 the stored half holds the conjugate of
            ! the mirrored mode, and the sign folds that into the index map for
            ! the other two axes.
            if (k1 >= 0) then
               j1 = k1; j2 = modulo(k2, n2); conjugateSign = 1
            else
               j1 = -k1; j2 = modulo(-k2, n2); conjugateSign = -1
            end if

            do a3 = 0, normalModes - 1
               k3 = mode_of_bin(a3, normalModes)
               j3 = modulo(conjugateSign*k3, n3)
               deconvolution = windowFactor1(j1)*windowFactor2(a2)*windowFactor3(a3)
               magnitudeSquared = real(spectrum(j1, j2, j3))**2 &
                                  + aimag(spectrum(j1, j2, j3))**2
               energy = energy + kernelCoefficients(a3)*deconvolution*deconvolution &
                        *magnitudeSquared
            end do
         end do
      end do
      energy = 0.5_dp*energy
      call system_clock(tick3)

      nullify (spectrum, spectrumBase)

      if (present(phaseTimes)) &
         phaseTimes = [real(tick1 - tick0, dp), real(tick2 - tick1, dp), &
                       0.0_dp, real(tick3 - tick2, dp)]/real(tickRate, dp)

   end function slab_energy

   !> Window deconvolution factors for one axis, indexed by FFT bin.  They are
   !> real and even in the mode index.
   subroutine window_factors(nModes, n, b, factors)

      !> Number of retained modes along this axis.
      integer, intent(in) :: nModes

      !> Fine grid extent along this axis.
      integer, intent(in) :: n

      !> Window shape parameter.
      real(dp), intent(in) :: b

      !> One factor per bin.
      real(dp), intent(out) :: factors(0:nModes - 1)

      integer :: a

      do a = 0, nModes - 1
         factors(a) = exp(b*(pi*real(mode_of_bin(a, nModes), dp)/real(n, dp))**2)
      end do

   end subroutine window_factors

   !> The same factors for the axis whose spectrum is stored as a Hermitian
   !> half.  Being even in the mode index, they can be tabulated by |k| alone.
   subroutine window_factors_half(nModes, n, b, factors)

      !> Number of retained modes along this axis.
      integer, intent(in) :: nModes

      !> Fine grid extent along this axis.
      integer, intent(in) :: n

      !> Window shape parameter.
      real(dp), intent(in) :: b

      !> One factor per non-negative mode index.
      real(dp), intent(out) :: factors(0:nModes/2)

      integer :: j

      do j = 0, nModes/2
         factors(j) = exp(b*(pi*real(j, dp)/real(n, dp))**2)
      end do

   end subroutine window_factors_half

   !> Fourier contribution to the per-atom potential and force of a monolayer.
   !>
   !> The kernel is a scalar per in-plane mode, so the structure is that of the
   !> fully periodic method: one adjoint transform for the structure factors,
   !> one forward transform of the weighted coefficients for the potential, and
   !> one forward transform per force component.
   !>
   !> The force normal to the layer comes out zero, exactly rather than by
   !> imposition: a coplanar charge distribution is symmetric under reflection
   !> in its own plane, so the potential is even in the normal coordinate.  In
   !> the code that shows up as a component whose wavevector is identically zero
   !> over the whole mode set, and its transform is skipped.
   subroutine long_range_potential_force_monolayer(positions, charges, nParticle, &
                                                   recVecs, area, params, energy, pot, &
                                                   force)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> In-plane crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(2, 3)

      !> Cell area.
      real(dp), intent(in) :: area

      !> Parameters of this evaluation.
      type(TEwaldParameters2d), intent(in) :: params

      !> Fourier contribution to the total energy.
      real(dp), intent(out) :: energy

      !> Fourier part of the potential at each charge.
      real(dp), intent(out) :: pot(:)

      !> Fourier part of the force on each charge, shape (3, nParticle).
      real(dp), intent(out), optional :: force(:, :)

      integer  :: nModes(2)
      real(dp), allocatable :: t1(:), t2(:)
      complex(dp), allocatable :: coefficients(:, :)   ! Shat, then weighted
      complex(dp), allocatable :: derivative(:, :)
      complex(dp), allocatable :: atSites(:)
      integer  :: i, a1, a2, k1, k2, iComponent
      real(dp) :: waveVector(3)
      real(dp) :: wavenumber
      real(dp) :: modeWeight
      real(dp) :: largestComponent(3)   ! largest |k_d| over the retained modes

      nModes = params%nModes(1:2)

      allocate (t1(nParticle), t2(nParticle))
      do i = 1, nParticle
         t1(i) = dot_product(recVecs(1, :), positions(i, :))
         t2(i) = dot_product(recVecs(2, :), positions(i, :))
      end do

      allocate (coefficients(0:nModes(1) - 1, 0:nModes(2) - 1))
      call adjoint_nfft_2d(nParticle, t1, t2, charges(1:nParticle), nModes(1), nModes(2), &
                           params%windowCutoff, coefficients)

      ! Weight the structure factors in place, taking the energy from the same
      ! pass so that no second array of this size is needed.
      energy = 0.0_dp
      largestComponent = 0.0_dp
      do a2 = 0, nModes(2) - 1
         k2 = mode_of_bin(a2, nModes(2))
         do a1 = 0, nModes(1) - 1
            k1 = mode_of_bin(a1, nModes(1))
            waveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                    + real(k2, dp)*recVecs(2, :))
            wavenumber = sqrt(sum(waveVector**2))
            modeWeight = regularised_kernel(wavenumber, 0.0_dp, area, params%alpha)
            energy = energy + modeWeight*(real(coefficients(a1, a2))**2 &
                                          + aimag(coefficients(a1, a2))**2)
            coefficients(a1, a2) = modeWeight*coefficients(a1, a2)
            largestComponent = max(largestComponent, abs(waveVector))
         end do
      end do
      energy = 0.5_dp*energy

      allocate (atSites(nParticle))
      call forward_nfft_2d(nParticle, t1, t2, nModes(1), nModes(2), params%windowCutoff, &
                           coefficients, atSites)
      pot(1:nParticle) = real(atSites(1:nParticle), dp)

      if (.not. present(force)) return

      force(:, 1:nParticle) = 0.0_dp
      allocate (derivative(0:nModes(1) - 1, 0:nModes(2) - 1))
      do iComponent = 1, 3
         ! A component with no wavevector anywhere in the mode set contributes
         ! nothing; for the usual geometry that is the normal component.
         if (largestComponent(iComponent) == 0.0_dp) cycle
         do a2 = 0, nModes(2) - 1
            k2 = mode_of_bin(a2, nModes(2))
            do a1 = 0, nModes(1) - 1
               k1 = mode_of_bin(a1, nModes(1))
               waveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                       + real(k2, dp)*recVecs(2, :))
               derivative(a1, a2) = cmplx(0.0_dp, waveVector(iComponent), dp) &
                                    *coefficients(a1, a2)
            end do
         end do
         call forward_nfft_2d(nParticle, t1, t2, nModes(1), nModes(2), &
                              params%windowCutoff, derivative, atSites)
         do i = 1, nParticle
            force(iComponent, i) = charges(i)*real(atSites(i), dp)
         end do
      end do

   end subroutine long_range_potential_force_monolayer

   !> Fourier contribution to the per-atom potential and force of a slab.
   !>
   !> The monolayer path with the open direction restored: the kernel is the
   !> periodic regularised one rather than its value at zero separation, so its
   !> coefficients depend on the mode along that direction and one small
   !> transform of the kernel is needed per in-plane mode, as in the energy.
   !> The wavevector the force differentiates to now has all three components.
   !>
   !> The regularisation parameters do not enter the derivative.  The period is
   !> built from the computed thickness so that every separation that occurs
   !> lies in the untouched region of the kernel, which makes the energy
   !> independent of the period and of the layer's centre; their derivatives are
   !> therefore zero and are not evaluated. Residual dependence is limited by
   !> the series truncation.
   !>
   !> The force is one order less smooth than the energy: the regularised kernel
   !> is only p-1 times differentiable across the joint, so its coefficients
   !> decay one power more slowly once differentiated, and the normal component
   !> pays most because its multiplier is largest at the edge of the band.  The
   !> shipped smoothness and mode-count floor are calibrated on the energy.
   !>
   !> Unlike a monolayer, a slab has no reflection symmetry, so its normal force
   !> is genuinely non-zero and all three components are transformed.
   subroutine long_range_potential_force_slab(positions, charges, nParticle, recVecs, &
                                              area, params, energy, pot, force)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> In-plane crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(2, 3)

      !> Cell area.
      real(dp), intent(in) :: area

      !> Parameters of this evaluation, including the period and centre the
      !> regularisation was sized for.
      type(TEwaldParameters2d), intent(in) :: params

      !> Fourier contribution to the total energy.
      real(dp), intent(out) :: energy

      !> Fourier part of the potential at each charge.
      real(dp), intent(out) :: pot(:)

      !> Fourier part of the force on each charge, shape (3, nParticle).
      real(dp), intent(out), optional :: force(:, :)

      real(dp), allocatable :: t1(:), t2(:), t3(:)
      real(dp), allocatable :: kernelCoefficients(:)
      complex(dp), allocatable :: coefficients(:, :, :)
      complex(dp), allocatable :: derivative(:, :, :)
      complex(dp), allocatable :: atSites(:)
      integer  :: i, a1, a2, a3, k1, k2, k3, iComponent
      integer  :: nModes(3)
      real(dp) :: period
      real(dp) :: inPlaneWaveVector(3)
      real(dp) :: fullWaveVector(3)      ! with the open-direction component added
      real(dp) :: wavenumber
      real(dp) :: coefficient
      real(dp) :: largestComponent(3)

      nModes = params%nModes
      period = params%period

      allocate (t1(nParticle), t2(nParticle), t3(nParticle))
      do i = 1, nParticle
         t1(i) = dot_product(recVecs(1, :), positions(i, :))
         t2(i) = dot_product(recVecs(2, :), positions(i, :))
         t3(i) = (positions(i, 3) - params%centreZ)/period
      end do

      allocate (coefficients(0:nModes(1) - 1, 0:nModes(2) - 1, 0:nModes(3) - 1))
      call adjoint_nfft_3d(nParticle, t1, t2, t3, charges(1:nParticle), nModes(1), &
                           nModes(2), nModes(3), params%windowCutoff, coefficients)

      ! Weight the structure factors in place, taking the energy from the same
      ! pass, as the monolayer path does.  The largest component of the
      ! wavevector is recorded so that a direction with no wavevector at all
      ! costs no transform below.
      allocate (kernelCoefficients(0:nModes(3) - 1))
      energy = 0.0_dp
      largestComponent = 0.0_dp

      do a2 = 0, nModes(2) - 1
         k2 = mode_of_bin(a2, nModes(2))
         do a1 = 0, nModes(1) - 1
            k1 = mode_of_bin(a1, nModes(1))
            inPlaneWaveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                           + real(k2, dp)*recVecs(2, :))
            wavenumber = sqrt(sum(inPlaneWaveVector**2))
            call kernel_fourier_coefficients(wavenumber, area, params%alpha, period, &
                                             params%layerWidth, params%smoothness, &
                                             nModes(3), kernelCoefficients)

            do a3 = 0, nModes(3) - 1
               k3 = mode_of_bin(a3, nModes(3))
               coefficient = kernelCoefficients(a3)
               energy = energy + coefficient*(real(coefficients(a1, a2, a3))**2 &
                                              + aimag(coefficients(a1, a2, a3))**2)
               coefficients(a1, a2, a3) = coefficient*coefficients(a1, a2, a3)

               fullWaveVector = inPlaneWaveVector
               fullWaveVector(3) = fullWaveVector(3) + 2.0_dp*pi*real(k3, dp)/period
               largestComponent = max(largestComponent, abs(fullWaveVector))
            end do
         end do
      end do
      energy = 0.5_dp*energy

      allocate (atSites(nParticle))
      call forward_nfft_3d(nParticle, t1, t2, t3, nModes(1), nModes(2), nModes(3), &
                           params%windowCutoff, coefficients, atSites)
      pot(1:nParticle) = real(atSites(1:nParticle), dp)

      if (.not. present(force)) return

      force(:, 1:nParticle) = 0.0_dp
      allocate (derivative(0:nModes(1) - 1, 0:nModes(2) - 1, 0:nModes(3) - 1))
      do iComponent = 1, 3
         if (largestComponent(iComponent) == 0.0_dp) cycle
         do a2 = 0, nModes(2) - 1
            k2 = mode_of_bin(a2, nModes(2))
            do a1 = 0, nModes(1) - 1
               k1 = mode_of_bin(a1, nModes(1))
               inPlaneWaveVector = 2.0_dp*pi*(real(k1, dp)*recVecs(1, :) &
                                              + real(k2, dp)*recVecs(2, :))
               do a3 = 0, nModes(3) - 1
                  k3 = mode_of_bin(a3, nModes(3))
                  fullWaveVector = inPlaneWaveVector
                  fullWaveVector(3) = fullWaveVector(3) + 2.0_dp*pi*real(k3, dp)/period
                  derivative(a1, a2, a3) = cmplx(0.0_dp, fullWaveVector(iComponent), dp) &
                                           *coefficients(a1, a2, a3)
               end do
            end do
         end do
         call forward_nfft_3d(nParticle, t1, t2, t3, nModes(1), nModes(2), nModes(3), &
                              params%windowCutoff, derivative, atSites)
         do i = 1, nParticle
            force(iComponent, i) = charges(i)*real(atSites(i), dp)
         end do
      end do

   end subroutine long_range_potential_force_slab

end module ewald_fft_2d_fourier

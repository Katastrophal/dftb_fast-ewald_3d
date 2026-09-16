module ewald_direct_2d
   !> Direct Ewald summation for a cell periodic in two directions and open in
   !> the third.
   !>
   !>   U = U_real + U_fourier + U_zeroMode + U_self
   !>
   !> The Fourier kernel depends on the pair separation along the open direction,
   !> so this branch cannot be factorised into squared structure factors. The
   !> zero in-plane Fourier coefficient is evaluated separately. Only
   !> charge-neutral configurations are accepted.
   use ewald_constants, only: dp, pi, sqrt_pi
   use ewald_geometry, only: TCell2d, cell_metrics_2d
   use ewald_truncation, only: direct_cutoffs
   use ewald_alpha_search, only: balanced_alpha, alphaSearchStatus
   use ewald_validation, only: check_configuration, minSeparationSquared, &
                               minWavenumberSquared
   use ewald_self, only: self_energy, self_potential
   use iso_fortran_env, only: error_unit
   implicit none

   private
   public :: ewald_energy, ewald_potential_force
   public :: real_space_energy, fourier_energy, zero_mode_energy
   public :: real_space_potential_force
   public :: balanced_splitting

   !> Accuracy used when the caller does not request one.
   real(dp), parameter :: defaultTolerance = 1.0e-12_dp

contains

   ! =====================================================================
   !  Public entry point
   ! =====================================================================

   !> Total electrostatic energy of the two-dimensionally periodic cell.
   !>
   !> Everything the method needs beyond the configuration itself is derived
   !> from a single requested accuracy: the splitting parameter is chosen to
   !> balance the two branches, and the two cutoffs follow from it.  Any of the
   !> three can be overridden, which is what the verification tests use to check
   !> that the energy does not depend on the splitting.
   function ewald_energy(positions, charges, nParticle, latVecs, tol, &
                         alpha_in, r_cut_in, k_cut_in) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.  The third
      !> component is the coordinate along the open direction and is not
      !> wrapped by anything.
      real(dp), intent(in) :: positions(:, :)

      !> Charges, in units where the Coulomb prefactor is one.
      real(dp), intent(in) :: charges(:)

      !> Number of charges in the cell.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.  Rows 1 and 2 span the periodic plane;
      !> row 3 is ignored.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy of the two truncations.  Defaults to 1e-12.
      real(dp), intent(in), optional :: tol

      !> Override for the splitting parameter.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the in-plane Fourier cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      !> Total energy of the cell.
      real(dp) :: energy

      type(TCell2d) :: cell
      real(dp) :: tolerance
      real(dp) :: alpha
      real(dp) :: r_cut, k_cut

      call check_configuration("ewald_direct_2d", positions, charges, nParticle)

      cell = cell_metrics_2d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_direct_2d: tol must lie in (0,1)"

      ! The search is skipped, not merely overwritten, when the caller pins the
      ! splitting: it is the expensive part of setting the method up, and a
      ! failed search would otherwise warn about a value that is discarded.
      if (present(alpha_in)) then
         alpha = alpha_in
      else
         alpha = balanced_splitting(cell, tolerance)
      end if
      if (alpha <= 0.0_dp) error stop "ewald_direct_2d: alpha must be positive"

      ! Deriving the cutoffs is skipped only if the caller pinned both of them.
      if (.not. (present(r_cut_in) .and. present(k_cut_in))) then
         call direct_cutoffs(alpha, tolerance, r_cut, k_cut)
      end if
      if (present(r_cut_in)) r_cut = r_cut_in
      if (present(k_cut_in)) k_cut = k_cut_in
      if (r_cut <= 0.0_dp .or. k_cut <= 0.0_dp) &
         error stop "ewald_direct_2d: cutoffs must be positive"

      energy = real_space_energy(positions, charges, nParticle, cell, alpha, r_cut) &
               + fourier_energy(positions, charges, nParticle, cell, alpha, k_cut) &
               + zero_mode_energy(positions, charges, nParticle, cell, alpha) &
               + self_energy(charges, nParticle, alpha)

   end function ewald_energy

   ! =====================================================================
   !  Real-space, Fourier and zero-mode branches
   ! =====================================================================

   !> Real-space branch: the screened pair interaction erfc(alpha*r)/r, summed
   !> over all pairs and over the in-plane periodic images inside the cutoff.
   !> Only two directions carry images, but the separation keeps all three
   !> components: the charges the images connect need not lie in one plane.
   !>
   !> Pairs are visited once each (j starts at i), which absorbs the factor of
   !> one half in front of the double sum.  A charge with its own images is
   !> visited once and therefore keeps that half explicitly.
   function real_space_energy(positions, charges, nParticle, cell, alpha, r_cut) &
      result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Geometry of the in-plane cell.
      type(TCell2d), intent(in) :: cell

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Cutoff radius beyond which pair terms are dropped.
      real(dp), intent(in) :: r_cut

      !> Real-space contribution to the total energy.
      real(dp) :: energy

      integer  :: nImages(2)              ! image shells along each direction
      integer  :: i, j, nx, ny
      integer  :: iDirection
      real(dp) :: dx, dy, dz
      real(dp) :: separationSquared, separation
      real(dp) :: cutoffSquared

      energy = 0.0_dp
      cutoffSquared = r_cut*r_cut

      ! The perpendicular width of the cell along direction i is 1/|b_i|, so
      ! the cutoff spans r_cut*|b_i| image shells; one extra is scanned so that
      ! no image is missed at the corners of the cell.
      do iDirection = 1, 2
         nImages(iDirection) = ceiling(r_cut*cell%recLengths(iDirection)) + 1
      end do

      !$omp parallel do default(shared) &
      !$omp private(j, nx, ny, dx, dy, dz, separationSquared, separation) &
      !$omp reduction(+:energy) schedule(dynamic)
      do i = 1, nParticle
         do j = i, nParticle
            do nx = -nImages(1), nImages(1)
               do ny = -nImages(2), nImages(2)
                  ! A charge does not interact with itself, but does interact
                  ! with all of its own images.
                  if (i == j .and. nx == 0 .and. ny == 0) cycle

                  dx = positions(i, 1) - positions(j, 1) &
                       + real(nx, dp)*cell%latVecs(1, 1) + real(ny, dp)*cell%latVecs(2, 1)
                  dy = positions(i, 2) - positions(j, 2) &
                       + real(nx, dp)*cell%latVecs(1, 2) + real(ny, dp)*cell%latVecs(2, 2)
                  dz = positions(i, 3) - positions(j, 3) &
                       + real(nx, dp)*cell%latVecs(1, 3) + real(ny, dp)*cell%latVecs(2, 3)

                  separationSquared = dx*dx + dy*dy + dz*dz

                  if (separationSquared < minSeparationSquared) &
                     error stop "ewald_direct_2d: coincident charges"
                  if (separationSquared > cutoffSquared) cycle
                  separation = sqrt(separationSquared)

                  if (i == j) then
                     energy = energy + 0.5_dp*charges(i)*charges(j) &
                              *erfc(alpha*separation)/separation
                  else
                     energy = energy + charges(i)*charges(j) &
                              *erfc(alpha*separation)/separation
                  end if
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end function real_space_energy

   !> Fourier branch: the smooth remainder, summed over in-plane modes with a
   !> non-zero wavevector.  Because the transform runs over two directions
   !> only, each mode carries a weight that still depends on the separation
   !> z_ij along the open direction:
   !>
   !>   Theta(k, z) = (pi/k) [ e^{ k z} erfc(k/(2 alpha) + alpha z)
   !>                        + e^{-k z} erfc(k/(2 alpha) - alpha z) ]
   !>
   !> The exponentials overflow and the error functions underflow long before
   !> their product does, so the pair is evaluated as erfc(x) = exp(-x^2)
   !> erfcx(x), which pulls a common Gaussian factor out in front and leaves
   !> only well-scaled quantities inside.
   !>
   !> The mode sum pairs k with -k, whose contributions are complex conjugates,
   !> so only the cosine of the in-plane phase survives.
   function fourier_energy(positions, charges, nParticle, cell, alpha, k_cut) &
      result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Geometry of the in-plane cell.
      type(TCell2d), intent(in) :: cell

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Cutoff wavenumber beyond which in-plane modes are dropped.
      real(dp), intent(in) :: k_cut

      !> Fourier contribution to the total energy.
      real(dp) :: energy

      integer  :: nModes(2)            ! mode indices to scan along each axis
      integer  :: mx, my
      integer  :: i, j
      integer  :: iDirection
      real(dp) :: kx, ky
      real(dp) :: wavenumberSquared, wavenumber
      real(dp) :: cutoffSquared
      real(dp) :: dx, dy
      real(dp) :: normalSeparation     ! z_ij, along the open direction
      real(dp) :: phase
      real(dp) :: gaussianFactor       ! the factor pulled out of both erfc terms
      real(dp) :: modeWeight           ! Theta(k, z_ij)

      energy = 0.0_dp
      cutoffSquared = k_cut*k_cut

      ! A mode index is m_i = k . a_i / (2*pi), so the cutoff admits
      ! k_cut*|a_i|/(2*pi) of them along direction i, plus one for the corners.
      do iDirection = 1, 2
         nModes(iDirection) = ceiling(k_cut*cell%latLengths(iDirection)/(2.0_dp*pi)) + 1
      end do

      !$omp parallel do default(shared) &
      !$omp private(my, kx, ky, wavenumberSquared, wavenumber, i, j) &
      !$omp private(dx, dy, normalSeparation) &
      !$omp private(phase, gaussianFactor, modeWeight) &
      !$omp reduction(+:energy) schedule(dynamic)
      do mx = -nModes(1), nModes(1)
         do my = -nModes(2), nModes(2)
            kx = (real(mx, dp)*cell%recVecs(1, 1) + real(my, dp)*cell%recVecs(2, 1))*2.0_dp*pi
            ky = (real(mx, dp)*cell%recVecs(1, 2) + real(my, dp)*cell%recVecs(2, 2))*2.0_dp*pi
            wavenumberSquared = kx*kx + ky*ky

            ! The zero mode is not dropped but handled separately, by
            ! zero_mode_energy; here it is the one mode that must be skipped.
            if (wavenumberSquared < minWavenumberSquared) cycle
            if (wavenumberSquared > cutoffSquared) cycle
            wavenumber = sqrt(wavenumberSquared)

            do i = 1, nParticle
               do j = 1, nParticle
                  dx = positions(i, 1) - positions(j, 1)
                  dy = positions(i, 2) - positions(j, 2)
                  normalSeparation = positions(i, 3) - positions(j, 3)
                  phase = kx*dx + ky*dy

                  gaussianFactor = exp(-wavenumberSquared/(4.0_dp*alpha**2) &
                                       - (alpha*normalSeparation)**2)
                  modeWeight = (pi/wavenumber)*gaussianFactor &
                               *(erfc_scaled(wavenumber/(2.0_dp*alpha) &
                                             + alpha*normalSeparation) &
                                 + erfc_scaled(wavenumber/(2.0_dp*alpha) &
                                               - alpha*normalSeparation))

                  energy = energy + charges(i)*charges(j)*cos(phase)*modeWeight
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! One factor of the area from the Fourier transform, one factor of two
      ! from the half in front of the double sum.
      energy = energy/(2.0_dp*cell%area)

   end function fourier_energy

   !> The mode at zero in-plane wavevector.
   !>
   !> In three dimensions this mode diverges and is dropped; here it stays
   !> finite and carries the interaction of each charge with the others smeared
   !> uniformly over sheets parallel to the layer.  For a strictly coplanar and
   !> neutral cell it vanishes: every separation along the open direction is
   !> zero, the weight becomes a constant, and the double sum collapses to the
   !> squared net charge.
   function zero_mode_energy(positions, charges, nParticle, cell, alpha) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Geometry of the cell; only the area is used here.
      type(TCell2d), intent(in) :: cell

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Contribution of the zero mode to the total energy.
      real(dp) :: energy

      integer  :: i, j
      real(dp) :: normalSeparation     ! z_ij, along the open direction
      real(dp) :: weight

      energy = 0.0_dp

      !$omp parallel do default(shared) private(j, normalSeparation, weight) &
      !$omp reduction(+:energy) schedule(dynamic)
      do i = 1, nParticle
         do j = 1, nParticle
            normalSeparation = positions(i, 3) - positions(j, 3)

            ! Zero-mode limit of the Fourier weight, written with sqrt(pi)
            ! divided into the second term so that the prefactor outside the
            ! sum is simply -pi/area.
            weight = normalSeparation*erf(alpha*normalSeparation) &
                     + exp(-alpha**2*normalSeparation**2)/(alpha*sqrt_pi)

            energy = energy + charges(i)*charges(j)*weight
         end do
      end do
      !$omp end parallel do

      energy = -(pi/cell%area)*energy

   end function zero_mode_energy

   ! =====================================================================
   !  Per-atom potential and force
   ! =====================================================================

   !> Electrostatic potential at each charge and the force acting on it, by the
   !> same explicit sums as the energy with one charge factor removed, so that
   !> U = (1/2) sum_i q_i phi(x_i) holds exactly.
   !>
   !> The force along the open direction is worth a note.  Differentiating the
   !> Fourier weight with respect to z_ij turns the sum of the two scaled error
   !> functions into their difference, the Gaussian pieces of the product rule
   !> cancelling exactly, so the derivative is as cheap as the weight itself.
   !> Both it and the zero mode's contribution are even in the separation, so
   !> both vanish term by term when every charge sits in one plane: the force
   !> normal to a monolayer is zero as an identity of the general formula, not
   !> as something imposed on it.
   subroutine ewald_potential_force(positions, charges, nParticle, latVecs, tol, &
                                    pot, force, energy, alpha_in, r_cut_in, k_cut_in)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy of the two truncations.  Defaults to 1e-12.
      real(dp), intent(in), optional :: tol

      !> Potential at each charge, pot(i) = phi(x_i).  Size at least nParticle.
      real(dp), intent(out) :: pot(:)

      !> Force on each charge, force(:, i) = -dU/dx_i.  Shape (3, nParticle).
      !> This is the physical force; a code that stores energy gradients wants
      !> its negative.  Optional: omitting it skips the force accumulation.
      real(dp), intent(out), optional :: force(:, :)

      !> Total energy, identical to what ewald_energy returns.
      real(dp), intent(out), optional :: energy

      !> Override for the splitting parameter.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the in-plane Fourier cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      type(TCell2d) :: cell
      logical  :: wantForce
      real(dp) :: tolerance
      real(dp) :: alpha, r_cut, k_cut
      real(dp) :: wavenumberCutoffSquared
      integer  :: nModes(2)
      integer  :: i, j, mx, my
      integer  :: iDirection
      real(dp) :: dx, dy
      real(dp) :: kx, ky
      real(dp) :: wavenumberSquared, wavenumber
      real(dp) :: normalSeparation        ! z_ij
      real(dp) :: phase
      real(dp) :: scaledArgumentPlus      ! k/(2 alpha) + alpha z
      real(dp) :: scaledArgumentMinus     ! k/(2 alpha) - alpha z
      real(dp) :: gaussianFactor
      real(dp) :: modeWeight              ! Theta(k, z)
      real(dp) :: modeWeightDerivative    ! dTheta/dz
      real(dp) :: modeForceFactor         ! scalar multiplying the wavevector
      real(dp) :: potentialAccumulator
      real(dp) :: forceX, forceY, forceZ

      call check_configuration("ewald_direct_2d", positions, charges, nParticle)
      wantForce = present(force)
      if (size(pot) < nParticle) &
         error stop "ewald_direct_2d: pot must be (N) and force (3,N)"
      if (wantForce) then
         if (size(force, 1) /= 3 .or. size(force, 2) < nParticle) &
            error stop "ewald_direct_2d: pot must be (N) and force (3,N)"
      end if

      cell = cell_metrics_2d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_direct_2d: tol must lie in (0,1)"

      if (present(alpha_in)) then
         alpha = alpha_in
      else
         alpha = balanced_splitting(cell, tolerance)
      end if
      if (alpha <= 0.0_dp) error stop "ewald_direct_2d: alpha must be positive"

      if (.not. (present(r_cut_in) .and. present(k_cut_in))) then
         call direct_cutoffs(alpha, tolerance, r_cut, k_cut)
      end if
      if (present(r_cut_in)) r_cut = r_cut_in
      if (present(k_cut_in)) k_cut = k_cut_in
      if (r_cut <= 0.0_dp .or. k_cut <= 0.0_dp) &
         error stop "ewald_direct_2d: cutoffs must be positive"

      wavenumberCutoffSquared = k_cut*k_cut
      do iDirection = 1, 2
         nModes(iDirection) = ceiling(k_cut*cell%latLengths(iDirection)/(2.0_dp*pi)) + 1
      end do

      ! --- real-space branch ------------------------------------------------
      ! A routine of its own because the fast method's small-cell fallback
      ! calls it too.  It leaves its contribution in pot and force, which the
      ! loop below picks up and adds the remaining branches to.
      call real_space_potential_force(positions, charges, nParticle, cell, &
                                      alpha, r_cut, pot, force)

      ! Each iteration owns one atom and writes only that atom's potential and
      ! force, so no reduction is needed.  Both remaining branches sit inside
      ! the loop, so the whole computation is parallel.
      !$omp parallel do default(shared) schedule(dynamic) &
      !$omp private(i, j, mx, my, dx, dy) &
      !$omp private(kx, ky) &
      !$omp private(wavenumberSquared, wavenumber) &
      !$omp private(normalSeparation, phase, scaledArgumentPlus, scaledArgumentMinus) &
      !$omp private(gaussianFactor, modeWeight, modeWeightDerivative) &
      !$omp private(modeForceFactor) &
      !$omp private(potentialAccumulator, forceX, forceY, forceZ)
      do i = 1, nParticle
         ! Resume from what the real-space branch left behind.
         potentialAccumulator = pot(i)
         forceX = 0.0_dp
         forceY = 0.0_dp
         forceZ = 0.0_dp
         if (wantForce) then
            forceX = force(1, i)
            forceY = force(2, i)
            forceZ = force(3, i)
         end if

         ! --- Fourier branch --------------------------------------------------
         do mx = -nModes(1), nModes(1)
            do my = -nModes(2), nModes(2)
               kx = (real(mx, dp)*cell%recVecs(1, 1) &
                     + real(my, dp)*cell%recVecs(2, 1))*2.0_dp*pi
               ky = (real(mx, dp)*cell%recVecs(1, 2) &
                     + real(my, dp)*cell%recVecs(2, 2))*2.0_dp*pi
               wavenumberSquared = kx*kx + ky*ky

               ! The zero mode is not dropped but handled by the block below.
               if (wavenumberSquared < minWavenumberSquared) cycle
               if (wavenumberSquared > wavenumberCutoffSquared) cycle
               wavenumber = sqrt(wavenumberSquared)

               do j = 1, nParticle
                  dx = positions(i, 1) - positions(j, 1)
                  dy = positions(i, 2) - positions(j, 2)
                  normalSeparation = positions(i, 3) - positions(j, 3)
                  phase = kx*dx + ky*dy

                  scaledArgumentPlus = wavenumber/(2.0_dp*alpha) + alpha*normalSeparation
                  scaledArgumentMinus = wavenumber/(2.0_dp*alpha) - alpha*normalSeparation
                  gaussianFactor = exp(-wavenumberSquared/(4.0_dp*alpha**2) &
                                       - (alpha*normalSeparation)**2)

                  modeWeight = (pi/wavenumber)*gaussianFactor &
                               *(erfc_scaled(scaledArgumentPlus) &
                                 + erfc_scaled(scaledArgumentMinus))

                  potentialAccumulator = potentialAccumulator &
                                         + charges(j)*cos(phase)*modeWeight/cell%area

                  if (wantForce) then
                     ! In-plane force: only the phase depends on the in-plane
                     ! coordinate, so differentiating it turns the cosine into
                     ! a sine and brings down the wavevector component.
                     modeForceFactor = charges(i)*charges(j)*sin(phase) &
                                       *modeWeight/cell%area
                     forceX = forceX + modeForceFactor*kx
                     forceY = forceY + modeForceFactor*ky

                     ! Normal force: here the weight itself depends on the
                     ! coordinate, and its derivative differs from the weight
                     ! only by the sign between the two scaled error functions.
                     modeWeightDerivative = pi*gaussianFactor &
                                            *(erfc_scaled(scaledArgumentPlus) &
                                              - erfc_scaled(scaledArgumentMinus))
                     forceZ = forceZ - charges(i)*charges(j)*cos(phase) &
                              *modeWeightDerivative/cell%area
                  end if
               end do
            end do
         end do

         ! --- zero mode --------------------------------------------------------
         do j = 1, nParticle
            normalSeparation = positions(i, 3) - positions(j, 3)
            potentialAccumulator = potentialAccumulator &
                                   - (2.0_dp*pi/cell%area)*charges(j) &
                                   *(normalSeparation*erf(alpha*normalSeparation) &
                                     + exp(-alpha**2*normalSeparation**2)/(alpha*sqrt_pi))
            if (wantForce) then
               forceZ = forceZ + (2.0_dp*pi/cell%area)*charges(i)*charges(j) &
                        *erf(alpha*normalSeparation)
            end if
         end do

         ! --- self term ---------------------------------------------------------
         pot(i) = potentialAccumulator + self_potential(charges(i), alpha)
         if (wantForce) then
            force(1, i) = forceX
            force(2, i) = forceY
            force(3, i) = forceZ
         end if
      end do
      !$omp end parallel do

      if (present(energy)) energy = 0.5_dp*dot_product(charges(1:nParticle), pot(1:nParticle))

   end subroutine ewald_potential_force

   !> Real-space branch of the per-atom potential and force: the terms
   !> real_space_energy sums, with one charge factor removed from the potential
   !> and accumulated per atom.  The loop over j therefore runs over all
   !> charges instead of starting at i, the pair sum no longer being halved.
   !>
   !> Public because the fast method calls it: when the plane is too small for
   !> the minimum-image convention its linked-cell path does not apply, and it
   !> falls back on this explicit image sum rather than carrying a second copy
   !> of it.
   subroutine real_space_potential_force(positions, charges, nParticle, cell, &
                                         alpha, r_cut, pot, force)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Geometry of the in-plane cell.
      type(TCell2d), intent(in) :: cell

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Cutoff radius beyond which pair terms are dropped.
      real(dp), intent(in) :: r_cut

      !> Real-space part of the potential at each charge.  Overwritten, not
      !> accumulated into.  Size at least nParticle.
      real(dp), intent(out) :: pot(:)

      !> Real-space part of the force on each charge, shape (3, nParticle).
      !> Overwritten, not accumulated into.  Optional: omitting it skips the
      !> force accumulation.
      real(dp), intent(out), optional :: force(:, :)

      logical  :: wantForce
      real(dp) :: cutoffSquared
      integer  :: nImages(2)
      integer  :: i, j, nx, ny
      integer  :: iDirection
      real(dp) :: dx, dy, dz
      real(dp) :: separationSquared, separation
      real(dp) :: complementaryError
      real(dp) :: pairForceFactor         ! scalar multiplying the separation vector
      real(dp) :: potentialAccumulator
      real(dp) :: forceX, forceY, forceZ

      wantForce = present(force)
      cutoffSquared = r_cut*r_cut

      ! Image shells, counted as in real_space_energy.
      do iDirection = 1, 2
         nImages(iDirection) = ceiling(r_cut*cell%recLengths(iDirection)) + 1
      end do

      ! Each iteration of the outer loop owns one atom and writes only that
      ! atom's potential and force, so no reduction is needed and the order in
      ! which the pair terms are summed does not depend on the thread count.
      !$omp parallel do default(shared) schedule(dynamic) &
      !$omp private(i, j, nx, ny, dx, dy, dz) &
      !$omp private(separationSquared, separation) &
      !$omp private(complementaryError, pairForceFactor) &
      !$omp private(potentialAccumulator, forceX, forceY, forceZ)
      do i = 1, nParticle
         potentialAccumulator = 0.0_dp
         forceX = 0.0_dp
         forceY = 0.0_dp
         forceZ = 0.0_dp

         do j = 1, nParticle
            do nx = -nImages(1), nImages(1)
               do ny = -nImages(2), nImages(2)
                  if (i == j .and. nx == 0 .and. ny == 0) cycle

                  dx = positions(i, 1) - positions(j, 1) &
                       + real(nx, dp)*cell%latVecs(1, 1) + real(ny, dp)*cell%latVecs(2, 1)
                  dy = positions(i, 2) - positions(j, 2) &
                       + real(nx, dp)*cell%latVecs(1, 2) + real(ny, dp)*cell%latVecs(2, 2)
                  dz = positions(i, 3) - positions(j, 3) &
                       + real(nx, dp)*cell%latVecs(1, 3) + real(ny, dp)*cell%latVecs(2, 3)

                  separationSquared = dx*dx + dy*dy + dz*dz
                  if (separationSquared < minSeparationSquared) &
                     error stop "ewald_direct_2d: coincident charges"
                  if (separationSquared > cutoffSquared) cycle
                  separation = sqrt(separationSquared)

                  complementaryError = erfc(alpha*separation)
                  potentialAccumulator = potentialAccumulator &
                                         + charges(j)*complementaryError/separation

                  ! Derivative of the screened pair term.  The separation
                  ! vector points from the image of j towards i, so a positive
                  ! factor pushes i away from j.
                  if (wantForce) then
                     pairForceFactor = charges(i)*charges(j) &
                                       *(2.0_dp*alpha/sqrt_pi &
                                         *exp(-alpha*alpha*separationSquared) &
                                         + complementaryError/separation) &
                                       /separationSquared
                     forceX = forceX + pairForceFactor*dx
                     forceY = forceY + pairForceFactor*dy
                     forceZ = forceZ + pairForceFactor*dz
                  end if
               end do
            end do
         end do

         pot(i) = potentialAccumulator
         if (wantForce) then
            force(1, i) = forceX
            force(2, i) = forceY
            force(3, i) = forceZ
         end if
      end do
      !$omp end parallel do

   end subroutine real_space_potential_force

   ! =====================================================================
   !  Splitting parameter
   ! =====================================================================

   !> Splitting parameter that makes the two branches converge at the same
   !> rate.  The bracketing and bisection are in ewald_alpha_search; what
   !> differs from the fully periodic case is that the shortest lattice and
   !> reciprocal vectors are taken over two directions and that a single
   !> Fourier mode contributes with the two-dimensional weight.
   function balanced_splitting(cell, tolerance) result(alpha)

      !> Geometry of the cell.
      type(TCell2d), intent(in) :: cell

      !> How closely the two branches have to be balanced.
      real(dp), intent(in) :: tolerance

      !> The splitting parameter.
      real(dp) :: alpha

      real(dp) :: shortestLatticeVector, shortestWaveVector
      integer  :: status

      shortestLatticeVector = minval(cell%latLengths)
      shortestWaveVector = 2.0_dp*pi*minval(cell%recLengths)

      call balanced_alpha(2, shortestWaveVector, shortestLatticeVector, cell%area, &
                 tolerance, alpha, status)

      if (status /= alphaSearchStatus%converged) then
         write (error_unit, '(A,I0)') &
            "ewald_direct_2d: could not balance the two branches, search status ", status
         ! Last resort: a fit of the balanced splitting against the cell size,
         ! with the area in place of a volume.  A heuristic, but the energy
         ! tolerates any positive value.
         alpha = exp(-0.310104_dp*log(cell%area) + 0.786382_dp)/2.0_dp
      end if

   end function balanced_splitting

end module ewald_direct_2d

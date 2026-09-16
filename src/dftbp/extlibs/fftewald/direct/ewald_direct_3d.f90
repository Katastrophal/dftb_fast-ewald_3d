module ewald_direct_3d
   !> Direct Ewald summation for a cell periodic in all three directions.
   !>
   !>   U = U_real + U_fourier + U_self
   !>
   !> Both branches are evaluated explicitly and scale quadratically. The
   !> Fourier sum omits the zero coefficient, corresponding to conducting
   !> boundary conditions at infinity. Only charge-neutral configurations are
   !> accepted.
   use ewald_constants, only: dp, pi, sqrt_pi
   use ewald_geometry, only: TCell3d, cell_metrics_3d
   use ewald_truncation, only: direct_cutoffs
   use ewald_alpha_search, only: balanced_alpha, alphaSearchStatus
   use ewald_validation, only: check_configuration, minSeparationSquared, &
                               minWavenumberSquared
   use ewald_self, only: self_energy, self_potential
   use iso_fortran_env, only: error_unit
   implicit none

   private
   public :: ewald_energy, ewald_potential_force
   public :: real_space_energy, fourier_energy
   public :: real_space_potential_force
   public :: balanced_splitting

   !> Accuracy used when the caller does not request one.
   real(dp), parameter :: defaultTolerance = 1.0e-12_dp

contains

   ! =====================================================================
   !  Public entry point
   ! =====================================================================

   !> Total electrostatic energy of the periodic cell.
   !>
   !> Everything beyond the configuration follows from the requested accuracy:
   !> the splitting parameter balances the two branches and the cutoffs follow
   !> from it.  Any of the three can be overridden, which is what the tests use
   !> to check that the energy does not depend on the splitting.
   function ewald_energy(positions, charges, nParticle, latVecs, tol, &
                         alpha_in, r_cut_in, k_cut_in) result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges, in units where the Coulomb prefactor is one.
      real(dp), intent(in) :: charges(:)

      !> Number of charges in the cell.
      integer, intent(in) :: nParticle

      !> Lattice vectors, one per row: latVecs(i, :) is a_i.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Requested accuracy of the two truncations.  Defaults to 1e-12.
      real(dp), intent(in), optional :: tol

      !> Override for the splitting parameter.  The energy is independent of
      !> it, so this only moves work between the two branches.
      real(dp), intent(in), optional :: alpha_in

      !> Override for the real-space cutoff radius.
      real(dp), intent(in), optional :: r_cut_in

      !> Override for the Fourier-space cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      !> Total energy of the cell.
      real(dp) :: energy

      type(TCell3d) :: cell
      real(dp) :: tolerance
      real(dp) :: alpha
      real(dp) :: r_cut, k_cut

      call check_configuration("ewald_direct_3d", positions, charges, nParticle)

      cell = cell_metrics_3d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_direct_3d: tol must lie in (0,1)"

      ! The search is skipped rather than overwritten when the caller pins the
      ! splitting: it is the expensive part of the set-up, and a failed search
      ! would warn about a value that is discarded.
      if (present(alpha_in)) then
         alpha = alpha_in
      else
         alpha = balanced_splitting(cell, tolerance)
      end if
      if (alpha <= 0.0_dp) error stop "ewald_direct_3d: alpha must be positive"

      ! Deriving the cutoffs is skipped only if the caller pinned both of them.
      if (.not. (present(r_cut_in) .and. present(k_cut_in))) then
         call direct_cutoffs(alpha, tolerance, r_cut, k_cut)
      end if
      if (present(r_cut_in)) r_cut = r_cut_in
      if (present(k_cut_in)) k_cut = k_cut_in
      if (r_cut <= 0.0_dp .or. k_cut <= 0.0_dp) &
         error stop "ewald_direct_3d: cutoffs must be positive"

      energy = real_space_energy(positions, charges, nParticle, cell, alpha, r_cut) &
               + fourier_energy(positions, charges, nParticle, cell, alpha, k_cut) &
               + self_energy(charges, nParticle, alpha)

   end function ewald_energy

   ! =====================================================================
   !  Real-space and Fourier branches
   ! =====================================================================

   !> Real-space branch: the screened pair interaction erfc(alpha*r)/r, summed
   !> over all pairs and over every periodic image inside the cutoff.
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

      !> Geometry of the cell.
      type(TCell3d), intent(in) :: cell

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Cutoff radius beyond which pair terms are dropped.
      real(dp), intent(in) :: r_cut

      !> Real-space contribution to the total energy.
      real(dp) :: energy

      integer  :: nImages(3)           ! image shells to scan per direction
      integer  :: i, j, nx, ny, nz
      real(dp) :: dx, dy, dz
      real(dp) :: separationSquared, separation
      real(dp) :: cutoffSquared
      integer  :: iDirection

      energy = 0.0_dp
      cutoffSquared = r_cut*r_cut

      ! The perpendicular width of the cell along direction i is 1/|b_i|, so
      ! the cutoff spans r_cut*|b_i| image shells; one extra is scanned so that
      ! no image is missed at the corners of the box.
      do iDirection = 1, 3
         nImages(iDirection) = ceiling(r_cut*cell%recLengths(iDirection)) + 1
      end do

      !$omp parallel do default(shared) &
      !$omp private(j, nx, ny, nz, dx, dy, dz, separationSquared, separation) &
      !$omp reduction(+:energy) schedule(dynamic)
      do i = 1, nParticle
         do j = i, nParticle
            do nx = -nImages(1), nImages(1)
               do ny = -nImages(2), nImages(2)
                  do nz = -nImages(3), nImages(3)
                     ! A charge does not interact with itself, but does
                     ! interact with all of its own images.
                     if (i == j .and. nx == 0 .and. ny == 0 .and. nz == 0) cycle

                     dx = positions(i, 1) - positions(j, 1) &
                          + real(nx, dp)*cell%latVecs(1, 1) &
                          + real(ny, dp)*cell%latVecs(2, 1) &
                          + real(nz, dp)*cell%latVecs(3, 1)
                     dy = positions(i, 2) - positions(j, 2) &
                          + real(nx, dp)*cell%latVecs(1, 2) &
                          + real(ny, dp)*cell%latVecs(2, 2) &
                          + real(nz, dp)*cell%latVecs(3, 2)
                     dz = positions(i, 3) - positions(j, 3) &
                          + real(nx, dp)*cell%latVecs(1, 3) &
                          + real(ny, dp)*cell%latVecs(2, 3) &
                          + real(nz, dp)*cell%latVecs(3, 3)

                     separationSquared = dx*dx + dy*dy + dz*dz

                     if (separationSquared < minSeparationSquared) &
                        error stop "ewald_direct_3d: coincident charges"
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
      end do
      !$omp end parallel do

   end function real_space_energy

   !> Fourier branch: the smooth remainder of the splitting, summed over
   !> reciprocal lattice modes inside the cutoff.  Each mode contributes the
   !> weighted squared magnitude of the structure factor
   !> S(k) = sum_j q_j exp(i k . x_j).
   !>
   !> Only half of the mode set is visited: S(-k) is the conjugate of S(k) and
   !> the weight is even, so the pair (k, -k) contributes twice what k does and
   !> the factor of one half in front of the sum cancels against it.
   !>
   !> The k = 0 mode is excluded, which is the conducting boundary condition.
   function fourier_energy(positions, charges, nParticle, cell, alpha, k_cut) &
      result(energy)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Charges.
      real(dp), intent(in) :: charges(:)

      !> Number of charges.
      integer, intent(in) :: nParticle

      !> Geometry of the cell.
      type(TCell3d), intent(in) :: cell

      !> Splitting parameter.
      real(dp), intent(in) :: alpha

      !> Cutoff wavenumber beyond which modes are dropped.
      real(dp), intent(in) :: k_cut

      !> Fourier-space contribution to the total energy.
      real(dp) :: energy

      integer  :: nModes(3)            ! mode indices to scan along each axis
      integer  :: mx, my, mz
      integer  :: j
      integer  :: iDirection
      real(dp) :: kx, ky, kz
      real(dp) :: wavenumberSquared
      real(dp) :: cutoffSquared
      real(dp) :: phase
      complex(dp) :: structureFactor

      energy = 0.0_dp
      cutoffSquared = k_cut*k_cut

      ! A mode index is m_i = k . a_i / (2*pi), so the cutoff admits
      ! k_cut*|a_i|/(2*pi) of them along direction i, plus one for the corners.
      do iDirection = 1, 3
         nModes(iDirection) = ceiling(k_cut*cell%latLengths(iDirection)/(2.0_dp*pi)) + 1
      end do

      !$omp parallel do default(shared) &
      !$omp private(my, mz, kx, ky, kz, wavenumberSquared, structureFactor, phase, j) &
      !$omp reduction(+:energy) schedule(dynamic)
      do mx = -nModes(1), nModes(1)
         do my = -nModes(2), nModes(2)
            do mz = 0, nModes(3)
               ! Keep one representative of each (k, -k) pair: the half space
               ! mz > 0, plus half of the plane mz = 0, plus half of the line
               ! mz = my = 0.
               if (mz == 0 .and. my < 0) cycle
               if (mz == 0 .and. my == 0 .and. mx <= 0) cycle

               kx = (real(mx, dp)*cell%recVecs(1, 1) + real(my, dp)*cell%recVecs(2, 1) &
                     + real(mz, dp)*cell%recVecs(3, 1))*2.0_dp*pi
               ky = (real(mx, dp)*cell%recVecs(1, 2) + real(my, dp)*cell%recVecs(2, 2) &
                     + real(mz, dp)*cell%recVecs(3, 2))*2.0_dp*pi
               kz = (real(mx, dp)*cell%recVecs(1, 3) + real(my, dp)*cell%recVecs(2, 3) &
                     + real(mz, dp)*cell%recVecs(3, 3))*2.0_dp*pi
               wavenumberSquared = kx*kx + ky*ky + kz*kz

               if (wavenumberSquared < minWavenumberSquared) cycle
               if (wavenumberSquared > cutoffSquared) cycle

               structureFactor = cmplx(0.0_dp, 0.0_dp, kind=dp)
               do j = 1, nParticle
                  phase = positions(j, 1)*kx + positions(j, 2)*ky + positions(j, 3)*kz
                  structureFactor = structureFactor &
                                    + charges(j)*exp(cmplx(0.0_dp, phase, kind=dp))
               end do

               energy = energy + exp(-wavenumberSquared/(4.0_dp*alpha**2)) &
                        /wavenumberSquared*abs(structureFactor)**2
            end do
         end do
      end do
      !$omp end parallel do

      energy = energy*4.0_dp*pi/cell%volume

   end function fourier_energy

   ! =====================================================================
   !  Per-atom potential and force
   ! =====================================================================

   !> Electrostatic potential at each charge and the force acting on it: the
   !> same explicit sums as the energy with one charge factor removed, so that
   !> U = (1/2) sum_i q_i phi(x_i) holds exactly.  O(N^2), like the energy.
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

      !> Override for the Fourier-space cutoff wavenumber.
      real(dp), intent(in), optional :: k_cut_in

      type(TCell3d) :: cell
      logical  :: wantForce
      real(dp) :: tolerance
      real(dp) :: alpha
      real(dp) :: r_cut, k_cut
      real(dp) :: wavenumberCutoffSquared
      integer  :: nModes(3)
      integer  :: i, j, mx, my, mz
      integer  :: iDirection
      real(dp) :: kx, ky, kz
      real(dp) :: wavenumberSquared
      real(dp) :: modeWeight           ! 4 pi exp(-k^2/4 alpha^2) / (V k^2)
      real(dp) :: modeForceFactor      ! scalar multiplying the wavevector
      real(dp) :: phase
      complex(dp) :: structureFactor
      complex(dp) :: localAmplitude    ! S(k) exp(-i k . x_i)

      call check_configuration("ewald_direct_3d", positions, charges, nParticle)
      wantForce = present(force)
      if (size(pot) < nParticle) &
         error stop "ewald_direct_3d: pot must be (N) and force (3,N)"
      if (wantForce) then
         if (size(force, 1) /= 3 .or. size(force, 2) < nParticle) &
            error stop "ewald_direct_3d: pot must be (N) and force (3,N)"
      end if

      cell = cell_metrics_3d(latVecs)

      tolerance = defaultTolerance
      if (present(tol)) tolerance = tol
      if (tolerance <= 0.0_dp .or. tolerance >= 1.0_dp) &
         error stop "ewald_direct_3d: tol must lie in (0,1)"

      if (present(alpha_in)) then
         alpha = alpha_in
      else
         alpha = balanced_splitting(cell, tolerance)
      end if
      if (alpha <= 0.0_dp) error stop "ewald_direct_3d: alpha must be positive"

      if (.not. (present(r_cut_in) .and. present(k_cut_in))) then
         call direct_cutoffs(alpha, tolerance, r_cut, k_cut)
      end if
      if (present(r_cut_in)) r_cut = r_cut_in
      if (present(k_cut_in)) k_cut = k_cut_in
      if (r_cut <= 0.0_dp .or. k_cut <= 0.0_dp) &
         error stop "ewald_direct_3d: cutoffs must be positive"

      wavenumberCutoffSquared = k_cut*k_cut

      ! --- real-space branch ------------------------------------------------
      ! A routine of its own because the fast method's small-cell fallback
      ! calls it too.
      call real_space_potential_force(positions, charges, nParticle, cell, &
                                      alpha, r_cut, pot, force)

      ! --- Fourier branch ---------------------------------------------------
      ! Unlike the energy, this needs the signed structure factor rather than
      ! its magnitude, so the full mode set is visited instead of half of it.
      do iDirection = 1, 3
         nModes(iDirection) = ceiling(k_cut*cell%latLengths(iDirection)/(2.0_dp*pi)) + 1
      end do

      do mx = -nModes(1), nModes(1)
         do my = -nModes(2), nModes(2)
            do mz = -nModes(3), nModes(3)
               kx = (real(mx, dp)*cell%recVecs(1, 1) + real(my, dp)*cell%recVecs(2, 1) &
                     + real(mz, dp)*cell%recVecs(3, 1))*2.0_dp*pi
               ky = (real(mx, dp)*cell%recVecs(1, 2) + real(my, dp)*cell%recVecs(2, 2) &
                     + real(mz, dp)*cell%recVecs(3, 2))*2.0_dp*pi
               kz = (real(mx, dp)*cell%recVecs(1, 3) + real(my, dp)*cell%recVecs(2, 3) &
                     + real(mz, dp)*cell%recVecs(3, 3))*2.0_dp*pi
               wavenumberSquared = kx*kx + ky*ky + kz*kz
               if (wavenumberSquared < minWavenumberSquared) cycle
               if (wavenumberSquared > wavenumberCutoffSquared) cycle

               modeWeight = 4.0_dp*pi/cell%volume &
                            *exp(-wavenumberSquared/(4.0_dp*alpha*alpha))/wavenumberSquared

               structureFactor = (0.0_dp, 0.0_dp)
               do j = 1, nParticle
                  phase = kx*positions(j, 1) + ky*positions(j, 2) + kz*positions(j, 3)
                  structureFactor = structureFactor &
                                    + charges(j)*cmplx(cos(phase), sin(phase), dp)
               end do

               do i = 1, nParticle
                  phase = kx*positions(i, 1) + ky*positions(i, 2) + kz*positions(i, 3)
                  localAmplitude = structureFactor*cmplx(cos(phase), -sin(phase), dp)
                  ! The real part of the local amplitude is the potential this
                  ! mode contributes; its imaginary part is what survives the
                  ! derivative of the phase, hence the force.
                  pot(i) = pot(i) + modeWeight*real(localAmplitude, dp)
                  if (wantForce) then
                     modeForceFactor = charges(i)*modeWeight*aimag(localAmplitude)
                     force(1, i) = force(1, i) - modeForceFactor*kx
                     force(2, i) = force(2, i) - modeForceFactor*ky
                     force(3, i) = force(3, i) - modeForceFactor*kz
                  end if
               end do
            end do
         end do
      end do

      ! --- self term --------------------------------------------------------
      ! It is independent of the positions, so it contributes no force.
      do i = 1, nParticle
         pot(i) = pot(i) + self_potential(charges(i), alpha)
      end do

      if (present(energy)) energy = 0.5_dp*dot_product(charges(1:nParticle), pot(1:nParticle))

   end subroutine ewald_potential_force

   !> Real-space branch of the per-atom potential and force: the terms
   !> real_space_energy sums, with one charge factor removed from the potential
   !> and accumulated per atom.  The loop over j therefore runs over all
   !> charges instead of starting at i, the pair sum no longer being halved.
   !>
   !> Public because the fast method calls it: when the cell is too small for
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

      !> Geometry of the cell.
      type(TCell3d), intent(in) :: cell

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
      integer  :: nImages(3)
      integer  :: i, j, nx, ny, nz
      integer  :: iDirection
      real(dp) :: dx, dy, dz
      real(dp) :: separationSquared, separation
      real(dp) :: pairForceFactor      ! scalar multiplying the separation vector
      real(dp) :: potentialAccumulator
      real(dp) :: forceX, forceY, forceZ

      wantForce = present(force)
      cutoffSquared = r_cut*r_cut

      ! Image shells, counted as in real_space_energy.
      do iDirection = 1, 3
         nImages(iDirection) = ceiling(r_cut*cell%recLengths(iDirection)) + 1
      end do

      ! Each iteration of the outer loop owns one atom and writes only that
      ! atom's potential and force, so no reduction is needed and the order in
      ! which the pair terms are summed does not depend on the thread count.
      !$omp parallel do default(shared) schedule(dynamic) &
      !$omp private(i, j, nx, ny, nz, dx, dy, dz, separationSquared) &
      !$omp private(separation, pairForceFactor) &
      !$omp private(potentialAccumulator, forceX, forceY, forceZ)
      do i = 1, nParticle
         potentialAccumulator = 0.0_dp
         forceX = 0.0_dp
         forceY = 0.0_dp
         forceZ = 0.0_dp

         do j = 1, nParticle
            do nx = -nImages(1), nImages(1)
               do ny = -nImages(2), nImages(2)
                  do nz = -nImages(3), nImages(3)
                     if (i == j .and. nx == 0 .and. ny == 0 .and. nz == 0) cycle

                     dx = positions(i, 1) - positions(j, 1) &
                          + real(nx, dp)*cell%latVecs(1, 1) &
                          + real(ny, dp)*cell%latVecs(2, 1) &
                          + real(nz, dp)*cell%latVecs(3, 1)
                     dy = positions(i, 2) - positions(j, 2) &
                          + real(nx, dp)*cell%latVecs(1, 2) &
                          + real(ny, dp)*cell%latVecs(2, 2) &
                          + real(nz, dp)*cell%latVecs(3, 2)
                     dz = positions(i, 3) - positions(j, 3) &
                          + real(nx, dp)*cell%latVecs(1, 3) &
                          + real(ny, dp)*cell%latVecs(2, 3) &
                          + real(nz, dp)*cell%latVecs(3, 3)

                     separationSquared = dx*dx + dy*dy + dz*dz
                     if (separationSquared < minSeparationSquared) &
                        error stop "ewald_direct_3d: coincident charges"
                     if (separationSquared > cutoffSquared) cycle
                     separation = sqrt(separationSquared)

                     potentialAccumulator = potentialAccumulator &
                                            + charges(j)*erfc(alpha*separation)/separation

                     ! Derivative of the screened pair term.  The separation
                     ! vector points from the image of j towards i, so a
                     ! positive factor pushes i away from j.
                     if (wantForce) then
                        pairForceFactor = charges(i)*charges(j) &
                                          *(2.0_dp*alpha/sqrt_pi &
                                            *exp(-alpha*alpha*separationSquared) &
                                            + erfc(alpha*separation)/separation) &
                                          /separationSquared
                        forceX = forceX + pairForceFactor*dx
                        forceY = forceY + pairForceFactor*dy
                        forceZ = forceZ + pairForceFactor*dz
                     end if
                  end do
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
   !> rate.  The bracketing and bisection are in ewald_alpha_search; what is
   !> supplied here is how much a single Fourier mode contributes in this
   !> geometry.
   function balanced_splitting(cell, tolerance) result(alpha)

      !> Geometry of the cell.
      type(TCell3d), intent(in) :: cell

      !> How closely the two branches have to be balanced.
      real(dp), intent(in) :: tolerance

      !> The splitting parameter.
      real(dp) :: alpha

      real(dp) :: shortestLatticeVector, shortestWaveVector
      integer  :: status

      shortestLatticeVector = minval(cell%latLengths)
      shortestWaveVector = 2.0_dp*pi*minval(cell%recLengths)

      call balanced_alpha(3, shortestWaveVector, shortestLatticeVector, cell%volume, &
                 tolerance, alpha, status)

      if (status /= alphaSearchStatus%converged) then
         write (error_unit, '(A,I0)') &
            "ewald_direct_3d: could not balance the two branches, search status ", status
         ! Last resort: a fit of the balanced splitting against the cell
         ! volume.  A heuristic, but the energy tolerates any positive value.
         alpha = exp(-0.310104_dp*log(cell%volume) + 0.786382_dp)/2.0_dp
      end if

   end function balanced_splitting

end module ewald_direct_3d

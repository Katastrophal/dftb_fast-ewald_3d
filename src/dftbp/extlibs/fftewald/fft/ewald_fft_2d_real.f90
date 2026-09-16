module ewald_fft_2d_real
   !> Real-space branch of the fast two-dimensional method.
   !>
   !> The linked cells tile the periodic plane. The separation along the open
   !> direction is retained without wrapping.
   !>
   !> When the plane is too small for the minimum-image convention the images
   !> are enumerated instead, by the direct reference's own sum, for the energy
   !> and for the potential and force alike.
   use ewald_constants, only: dp, sqrt_pi
   use ewald_geometry, only: TCell2d
   use ewald_cell_list, only: TCellList2d, build_cell_list_2d, box_counts_2d, &
                              box_index, minBoxesPerAxis
   use ewald_validation, only: minSeparationSquared
   use ewald_direct_2d, only: direct_real_space_energy => real_space_energy, &
                              direct_real_space_potential_force => real_space_potential_force
   implicit none

   private
   public :: short_range_energy, short_range_energy_cell_list
   public :: short_range_potential_force

contains

   !> Real-space contribution to the energy, by whichever path the cell allows.
   function short_range_energy(positions, charges, nParticle, cell, alpha, r_cut) &
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

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Real-space contribution to the total energy.
      real(dp) :: energy

      integer :: nBoxes(2)

      nBoxes = box_counts_2d(cell%recLengths, r_cut)

      if (minval(nBoxes) >= minBoxesPerAxis) then
         energy = short_range_energy_cell_list(positions, charges, nParticle, cell, &
                                               alpha, r_cut, nBoxes)
      else
         ! Too small for a unique minimum image: enumerate the in-plane images
         ! instead, using the direct reference's own real-space sum.
         energy = direct_real_space_energy(positions, charges, nParticle, cell, &
                                           alpha, r_cut)
      end if

   end function short_range_energy

   !> Real-space contribution to the energy through the linked-cell
   !> decomposition.
   function short_range_energy_cell_list(positions, charges, nParticle, cell, &
                                         alpha, r_cut, nBoxes) result(energy)

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

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Number of boxes along each in-plane axis.
      integer, intent(in) :: nBoxes(2)

      !> Real-space contribution to the total energy.
      real(dp) :: energy

      type(TCellList2d) :: list
      real(dp) :: a1(3), a2(3)         ! in-plane lattice vectors
      real(dp) :: cutoffSquared
      integer  :: i, j
      integer  :: ix, iy               ! box of charge i
      integer  :: dx, dy               ! offsets to a neighbouring box
      integer  :: jx, jy               ! box being scanned
      real(dp) :: ds1, ds2             ! in-plane fractional separation, nearest image
      real(dp) :: separationVector(3)
      real(dp) :: separationSquared, separation

      a1 = cell%latVecs(1, :)
      a2 = cell%latVecs(2, :)
      cutoffSquared = r_cut*r_cut

      call build_cell_list_2d(positions, nParticle, cell%recVecs, nBoxes, list)

      energy = 0.0_dp

      ! Each pair is met twice, once from each end, which the factor of one
      ! half accounts for.
      !$omp parallel do default(shared) &
      !$omp private(ix, iy, dx, dy, jx, jy, j, ds1, ds2) &
      !$omp private(separationVector, separationSquared, separation) &
      !$omp reduction(+:energy) schedule(guided)
      do i = 1, nParticle
         ix = box_index(list%frac1(i), nBoxes(1))
         iy = box_index(list%frac2(i), nBoxes(2))

         do dy = -1, 1
            jy = modulo(iy + dy, nBoxes(2))
            do dx = -1, 1
               jx = modulo(ix + dx, nBoxes(1))

               j = list%head(jx, jy)
               do while (j /= 0)
                  if (j /= i) then
                     ds1 = list%frac1(i) - list%frac1(j); ds1 = ds1 - nint(ds1)
                     ds2 = list%frac2(i) - list%frac2(j); ds2 = ds2 - nint(ds2)
                     separationVector = ds1*a1 + ds2*a2
                     ! The open direction carries no image, so its component is
                     ! the plain difference of the two coordinates.
                     separationVector(3) = separationVector(3) &
                                           + (positions(i, 3) - positions(j, 3))
                     separationSquared = sum(separationVector**2)

                     if (separationSquared <= cutoffSquared .and. &
                         separationSquared > minSeparationSquared) then
                        separation = sqrt(separationSquared)
                        energy = energy + 0.5_dp*charges(i)*charges(j) &
                                 *erfc(alpha*separation)/separation
                     end if
                  end if
                  j = list%next(j)
               end do

            end do
         end do
      end do
      !$omp end parallel do

   end function short_range_energy_cell_list

   !> Real-space contribution to the per-atom potential and force, by whichever
   !> path the cell allows.
   subroutine short_range_potential_force(positions, charges, nParticle, cell, &
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

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Real-space part of the potential at each charge.
      real(dp), intent(out) :: pot(:)

      !> Real-space part of the force on each charge, shape (3, nParticle).
      real(dp), intent(out) :: force(:, :)

      integer :: nBoxes(2)

      pot(1:nParticle) = 0.0_dp
      force(:, 1:nParticle) = 0.0_dp

      nBoxes = box_counts_2d(cell%recLengths, r_cut)

      if (minval(nBoxes) >= minBoxesPerAxis) then
         call potential_force_cell_list(positions, charges, nParticle, cell, &
                                        alpha, r_cut, nBoxes, pot, force)
      else
         ! The plane is too small for the minimum image to be unique, so the
         ! images have to be enumerated.  This is the direct reference's own
         ! per-atom real-space sum, reused rather than reimplemented.
         call direct_real_space_potential_force(positions, charges, nParticle, cell, &
                                                alpha, r_cut, pot, force)
      end if

   end subroutine short_range_potential_force

   !> Per-atom potential and force through the linked-cell decomposition.
   !>
   !> Mirrors short_range_energy_cell_list term for term, including its
   !> coincidence filter, so that the potential sums to exactly the energy that
   !> path returns.  Each iteration writes only its own atom, so no reduction
   !> is needed.
   subroutine potential_force_cell_list(positions, charges, nParticle, cell, &
                                        alpha, r_cut, nBoxes, pot, force)

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

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Number of boxes along each in-plane axis.
      integer, intent(in) :: nBoxes(2)

      !> Real-space part of the potential at each charge.
      real(dp), intent(out) :: pot(:)

      !> Real-space part of the force on each charge.
      real(dp), intent(out) :: force(:, :)

      type(TCellList2d) :: list
      real(dp) :: a1(3), a2(3)
      real(dp) :: cutoffSquared
      integer  :: i, j, ix, iy, dx, dy, jx, jy
      real(dp) :: ds1, ds2
      real(dp) :: separationVector(3)
      real(dp) :: separationSquared, separation
      real(dp) :: complementaryError
      real(dp) :: pairForceFactor
      real(dp) :: potentialAccumulator
      real(dp) :: forceAccumulator(3)

      a1 = cell%latVecs(1, :)
      a2 = cell%latVecs(2, :)
      cutoffSquared = r_cut*r_cut

      call build_cell_list_2d(positions, nParticle, cell%recVecs, nBoxes, list)

      !$omp parallel do default(shared) &
      !$omp private(i, ix, iy, dx, dy, jx, jy, j, ds1, ds2) &
      !$omp private(separationVector, separationSquared, separation) &
      !$omp private(complementaryError, pairForceFactor) &
      !$omp private(potentialAccumulator, forceAccumulator) schedule(guided)
      do i = 1, nParticle
         ix = box_index(list%frac1(i), nBoxes(1))
         iy = box_index(list%frac2(i), nBoxes(2))
         potentialAccumulator = 0.0_dp
         forceAccumulator = 0.0_dp

         do dy = -1, 1
            jy = modulo(iy + dy, nBoxes(2))
            do dx = -1, 1
               jx = modulo(ix + dx, nBoxes(1))

               j = list%head(jx, jy)
               do while (j /= 0)
                  if (j /= i) then
                     ds1 = list%frac1(i) - list%frac1(j); ds1 = ds1 - nint(ds1)
                     ds2 = list%frac2(i) - list%frac2(j); ds2 = ds2 - nint(ds2)
                     separationVector = ds1*a1 + ds2*a2
                     separationVector(3) = separationVector(3) &
                                           + (positions(i, 3) - positions(j, 3))
                     separationSquared = sum(separationVector**2)

                     if (separationSquared <= cutoffSquared .and. &
                         separationSquared > minSeparationSquared) then
                        separation = sqrt(separationSquared)
                        complementaryError = erfc(alpha*separation)
                        potentialAccumulator = potentialAccumulator &
                                               + charges(j)*complementaryError/separation
                        pairForceFactor = charges(i)*charges(j) &
                                          *(2.0_dp*alpha/sqrt_pi &
                                            *exp(-alpha*alpha*separationSquared) &
                                            + complementaryError/separation) &
                                          /separationSquared
                        forceAccumulator = forceAccumulator &
                                           + pairForceFactor*separationVector
                     end if
                  end if
                  j = list%next(j)
               end do

            end do
         end do

         pot(i) = potentialAccumulator
         force(:, i) = forceAccumulator
      end do
      !$omp end parallel do

   end subroutine potential_force_cell_list

end module ewald_fft_2d_real

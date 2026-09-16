module ewald_fft_3d_real
   !> Real-space branch of the fast three-dimensional method.
   !>
   !> The screened pair term is the one the direct reference sums; what changes
   !> is how the pairs are found.  Because the cutoff is tied to the mean
   !> spacing between charges, every partner within reach lies in a fixed
   !> neighbourhood of boxes, and a linked-cell decomposition finds them in
   !> O(N).  The branch is therefore exact up to the truncation at the cutoff,
   !> which the parameter recipe budgets for.
   !>
   !> The linked-cell path assumes the minimum-image convention, which requires
   !> at least three boxes along every axis.  A cell too small for that falls
   !> back on the direct reference's own explicit image sum; neither is
   !> reimplemented here.
   use ewald_constants, only: dp, sqrt_pi
   use ewald_geometry, only: TCell3d
   use ewald_cell_list, only: TCellList3d, build_cell_list_3d, box_counts_3d, &
                              box_index, minBoxesPerAxis
   use ewald_validation, only: minSeparationSquared
   use ewald_direct_3d, only: direct_real_space_potential_force => real_space_potential_force
   implicit none

   private
   public :: short_range_potential_force

contains

   !> Real-space contribution to the per-atom potential and force, by whichever
   !> path the cell allows.  The pair potential is differentiated for the force;
   !> the separation vector points from the image of the partner towards the
   !> charge, so a positive factor pushes the two apart.
   subroutine short_range_potential_force(positions, charges, nParticle, cell, &
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

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Real-space part of the potential at each charge.
      real(dp), intent(out) :: pot(:)

      !> Real-space part of the force on each charge, shape (3, nParticle).
      real(dp), intent(out) :: force(:, :)

      integer :: nBoxes(3)

      pot(1:nParticle) = 0.0_dp
      force(:, 1:nParticle) = 0.0_dp

      nBoxes = box_counts_3d(cell%recLengths, r_cut)

      if (minval(nBoxes) >= minBoxesPerAxis) then
         call potential_force_cell_list(positions, charges, nParticle, cell, &
                                        alpha, r_cut, nBoxes, pot, force)
      else
         ! The cell is too small for the minimum image to be unique, so the
         ! images have to be enumerated.  This is the direct reference's own
         ! per-atom real-space sum, reused rather than reimplemented.
         call direct_real_space_potential_force(positions, charges, nParticle, cell, &
                                                alpha, r_cut, pot, force)
      end if

   end subroutine short_range_potential_force

   !> Per-atom potential and force through the linked-cell decomposition.
   !> Each iteration writes only its own atom, so no reduction is needed.
   subroutine potential_force_cell_list(positions, charges, nParticle, cell, &
                                        alpha, r_cut, nBoxes, pot, force)

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

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Number of boxes along each axis.
      integer, intent(in) :: nBoxes(3)

      !> Real-space part of the potential at each charge.
      real(dp), intent(out) :: pot(:)

      !> Real-space part of the force on each charge.
      real(dp), intent(out) :: force(:, :)

      type(TCellList3d) :: list
      real(dp) :: a1(3), a2(3), a3(3)
      real(dp) :: cutoffSquared
      integer  :: i, j
      integer  :: ix, iy, iz, dx, dy, dz, jx, jy, jz
      real(dp) :: ds1, ds2, ds3
      real(dp) :: separationVector(3)
      real(dp) :: separationSquared, separation
      real(dp) :: complementaryError
      real(dp) :: pairForceFactor            ! scalar multiplying the separation vector
      real(dp) :: potentialAccumulator
      real(dp) :: forceAccumulator(3)

      a1 = cell%latVecs(1, :)
      a2 = cell%latVecs(2, :)
      a3 = cell%latVecs(3, :)
      cutoffSquared = r_cut*r_cut

      call build_cell_list_3d(positions, nParticle, cell%recVecs, nBoxes, list)

      !$omp parallel do default(shared) &
      !$omp private(i, ix, iy, iz, dx, dy, dz, jx, jy, jz, j, ds1, ds2, ds3) &
      !$omp private(separationVector, separationSquared, separation) &
      !$omp private(complementaryError, pairForceFactor) &
      !$omp private(potentialAccumulator, forceAccumulator) schedule(guided)
      do i = 1, nParticle
         ix = box_index(list%frac1(i), nBoxes(1))
         iy = box_index(list%frac2(i), nBoxes(2))
         iz = box_index(list%frac3(i), nBoxes(3))
         potentialAccumulator = 0.0_dp
         forceAccumulator = 0.0_dp

         do dz = -1, 1
            jz = modulo(iz + dz, nBoxes(3))
            do dy = -1, 1
               jy = modulo(iy + dy, nBoxes(2))
               do dx = -1, 1
                  jx = modulo(ix + dx, nBoxes(1))

                  j = list%head(jx, jy, jz)
                  do while (j /= 0)
                     if (j /= i) then
                        ds1 = list%frac1(i) - list%frac1(j); ds1 = ds1 - nint(ds1)
                        ds2 = list%frac2(i) - list%frac2(j); ds2 = ds2 - nint(ds2)
                        ds3 = list%frac3(i) - list%frac3(j); ds3 = ds3 - nint(ds3)
                        separationVector = ds1*a1 + ds2*a2 + ds3*a3
                        separationSquared = sum(separationVector**2)

                        if (separationSquared < minSeparationSquared) &
                           error stop "ewald_fft_3d: coincident charges"
                        if (separationSquared <= cutoffSquared) then
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
         end do

         pot(i) = potentialAccumulator
         force(:, i) = forceAccumulator
      end do
      !$omp end parallel do

   end subroutine potential_force_cell_list

end module ewald_fft_3d_real

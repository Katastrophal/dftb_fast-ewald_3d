module ewald_geometry
   !> Unit-cell geometry shared by every Ewald method in this library: cell size,
   !> dual basis and the lengths of both bases, one derived type per geometry.
   !>
   !> Conventions.  Lattice vectors are stored by rows, latVecs(i, :) being a_i;
   !> this is the transpose of the column convention some atomistic codes use.
   !> The stored dual basis is the crystallographic one, b_i . a_j = delta_ij,
   !> and carries no factor of 2*pi, so that a fractional coordinate is simply
   !> s_i = b_i . x.  A physical wavevector is formed where it is needed as
   !> k = 2*pi*(m1*b1 + m2*b2 + ...).
   use ewald_constants, only: dp
   implicit none

   private
   public :: cross_product
   public :: TCell3d, cell_metrics_3d

   !> Below these the cell counts as collapsed and its dual basis is meaningless.
   real(dp), parameter :: min_volume = 1.0e-18_dp

   !> Geometry of a cell that is periodic in all three directions.
   type :: TCell3d

      !> Real-space lattice vectors, one per row.
      real(dp) :: latVecs(3, 3) = 0.0_dp

      !> Crystallographic dual basis, one per row (no factor of 2*pi).
      real(dp) :: recVecs(3, 3) = 0.0_dp

      !> Cell volume |a1 . (a2 x a3)|.
      real(dp) :: volume = 0.0_dp

      !> Lengths |a_i|.  These, and not the perpendicular cell widths, set how
      !> far a mode index has to run to cover a given cutoff, a mode index
      !> being m_i = k . a_i / (2*pi).  The two coincide only for an orthogonal
      !> cell.
      real(dp) :: latLengths(3) = 0.0_dp

      !> Lengths |b_i|.  The perpendicular width of the cell along direction i
      !> is 1/|b_i|, which limits how wide a cutoff sphere may be before the
      !> minimum-image convention breaks down.
      real(dp) :: recLengths(3) = 0.0_dp

   end type TCell3d

contains

   !> Cross product of two vectors in three dimensions.
   pure function cross_product(a, b) result(c)

      !> Factors of the product.
      real(dp), intent(in) :: a(3), b(3)

      !> The product a x b.
      real(dp) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

   end function cross_product

   !> Derive the geometry of a three-dimensionally periodic cell from its
   !> lattice vectors.
   function cell_metrics_3d(latVecs) result(cell)

      !> Real-space lattice vectors, one per row.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Volume, dual basis and the lengths of both bases.
      type(TCell3d) :: cell

      cell%latVecs = latVecs

      ! The sign of the triple product only records the handedness of the basis.
      cell%volume = abs(dot_product(latVecs(1, :), &
                                    cross_product(latVecs(2, :), latVecs(3, :))))
      if (cell%volume <= min_volume) error stop "ewald_geometry: cell volume too small"

      ! b_i is perpendicular to the two lattice vectors other than a_i,
      ! normalised so that b_i . a_i = 1.
      cell%recVecs(1, :) = cross_product(latVecs(2, :), latVecs(3, :))/cell%volume
      cell%recVecs(2, :) = cross_product(latVecs(3, :), latVecs(1, :))/cell%volume
      cell%recVecs(3, :) = cross_product(latVecs(1, :), latVecs(2, :))/cell%volume

      cell%latLengths = [sqrt(sum(latVecs(1, :)**2)), &
                         sqrt(sum(latVecs(2, :)**2)), &
                         sqrt(sum(latVecs(3, :)**2))]
      cell%recLengths = [sqrt(sum(cell%recVecs(1, :)**2)), &
                         sqrt(sum(cell%recVecs(2, :)**2)), &
                         sqrt(sum(cell%recVecs(3, :)**2))]

   end function cell_metrics_3d

end module ewald_geometry

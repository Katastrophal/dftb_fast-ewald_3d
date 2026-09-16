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
   public :: TCell3d, TCell2d, cell_metrics_3d, cell_metrics_2d

   !> Below these the cell counts as collapsed and its dual basis is meaningless.
   real(dp), parameter :: min_volume = 1.0e-18_dp
   real(dp), parameter :: min_area = 1.0e-18_dp

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

   !> Geometry of a cell that is periodic in two directions and open in the
   !> third.  Rows 1 and 2 of latVecs span the periodic plane; row 3 is ignored
   !> and the extent along the normal is read off the particle positions.
   type :: TCell2d

      !> Real-space lattice vectors, one per row.  Only rows 1 and 2 are used.
      real(dp) :: latVecs(3, 3) = 0.0_dp

      !> In-plane crystallographic dual basis, one per row (no factor of 2*pi).
      real(dp) :: recVecs(2, 3) = 0.0_dp

      !> Area |a1 x a2| of the two-dimensional cell.
      real(dp) :: area = 0.0_dp

      !> The un-normalised plane normal a1 x a2, whose length is the area.
      real(dp) :: normal(3) = 0.0_dp

      !> Lengths |a_1|, |a_2| of the in-plane lattice vectors.
      real(dp) :: latLengths(2) = 0.0_dp

      !> Lengths |b_1|, |b_2| of the in-plane dual basis.
      real(dp) :: recLengths(2) = 0.0_dp

   end type TCell2d

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

   !> Derive the geometry of a two-dimensionally periodic cell.  Only the first
   !> two rows of latVecs are read.
   function cell_metrics_2d(latVecs) result(cell)

      !> Real-space lattice vectors, one per row; row 3 is ignored.
      real(dp), intent(in) :: latVecs(3, 3)

      !> Area, plane normal, dual basis and the lengths of both bases.
      type(TCell2d) :: cell

      real(dp) :: a1(3), a2(3)

      a1 = latVecs(1, :)
      a2 = latVecs(2, :)

      cell%latVecs = latVecs
      cell%normal = cross_product(a1, a2)
      cell%area = sqrt(sum(cell%normal**2))
      if (cell%area <= min_area) error stop "ewald_geometry: cell area too small"

      ! Crossing a2 (respectively a1) with the normal rotates it within the
      ! plane; dividing by the squared area normalises it to b_i . a_j = delta_ij.
      cell%recVecs(1, :) = cross_product(a2, cell%normal)/cell%area**2
      cell%recVecs(2, :) = cross_product(cell%normal, a1)/cell%area**2

      cell%latLengths = [sqrt(sum(a1**2)), sqrt(sum(a2**2))]
      cell%recLengths = [sqrt(sum(cell%recVecs(1, :)**2)), &
                         sqrt(sum(cell%recVecs(2, :)**2))]

   end function cell_metrics_2d

end module ewald_geometry

module ewald_cell_list
   !> Linked-cell decomposition of the simulation cell, used by the real-space
   !> branch of both fast methods.
   !>
   !> The cell is cut into boxes at least as wide as the cutoff, so every
   !> partner within reach of a charge lies in its own box or in one of the
   !> neighbouring ones.  That neighbourhood does not grow with N, which is what
   !> turns the O(N^2) pair sum into an O(N) one.  Storage is the classical
   !> head-and-next pair: head(cell) is one charge in that box, next(i) the
   !> charge after i in the same box, zero terminating the chain.
   !>
   !> The construction assumes the minimum-image convention, which holds
   !> precisely when at least three boxes fit along every periodic direction.
   !> Callers must check the box counts before building a list and fall back on
   !> an explicit sum over periodic images if the cell is smaller than that.
   use ewald_constants, only: dp
   implicit none

   private
   public :: TCellList3d, TCellList2d, build_cell_list_3d, build_cell_list_2d
   public :: box_index, box_counts_3d, box_counts_2d
   public :: minBoxesPerAxis

   !> Fewest boxes per periodic direction the decomposition needs for the
   !> minimum-image convention to hold: the charge's own box and one on each
   !> side.  A caller whose box counts fall below this must not build a list.
   integer, parameter :: minBoxesPerAxis = 3

   !> Linked-cell decomposition of a three-dimensionally periodic cell.
   type :: TCellList3d

      !> Number of boxes along each periodic direction.
      integer :: nBoxes(3) = 0

      !> Fractional coordinates wrapped into [0, 1), one array per lattice
      !> direction: the pair loop reads them one direction at a time.
      real(dp), allocatable :: frac1(:), frac2(:), frac3(:)

      !> head(ix, iy, iz) is one charge in that box, or zero if it is empty.
      integer, allocatable :: head(:, :, :)

      !> next(i) is the charge after i in the same box, or zero.
      integer, allocatable :: next(:)

   end type TCellList3d

   !> Linked-cell decomposition of a two-dimensionally periodic cell.  The
   !> boxes tile the periodic plane only: the open direction has no periodicity
   !> to exploit and its extent is small beside the in-plane cell.
   type :: TCellList2d

      !> Number of boxes along each in-plane direction.
      integer :: nBoxes(2) = 0

      !> In-plane fractional coordinates wrapped into [0, 1).
      real(dp), allocatable :: frac1(:), frac2(:)

      !> head(ix, iy) is one charge in that box, or zero if it is empty.
      integer, allocatable :: head(:, :)

      !> next(i) is the charge after i in the same box, or zero.
      integer, allocatable :: next(:)

   end type TCellList2d

contains

   !> Largest number of boxes that still leaves every box at least as wide as
   !> the cutoff, along each of three periodic directions.
   !>
   !> The perpendicular width of the cell along direction i is 1/|b_i|, so
   !> n boxes of width at least r_cut require n <= 1/(r_cut*|b_i|).  A result
   !> below three means the cell is too small for the minimum-image convention
   !> at this cutoff and the caller must not build a list.
   pure function box_counts_3d(recLengths, r_cut) result(nBoxes)

      !> Lengths |b_i| of the three dual-basis vectors.
      real(dp), intent(in) :: recLengths(3)

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Box counts along the three directions.
      integer :: nBoxes(3)

      integer :: iDirection

      do iDirection = 1, 3
         nBoxes(iDirection) = floor(1.0_dp/(r_cut*recLengths(iDirection)))
      end do

   end function box_counts_3d

   !> Two-dimensional counterpart of box_counts_3d.
   pure function box_counts_2d(recLengths, r_cut) result(nBoxes)

      !> Lengths |b_i| of the two in-plane dual-basis vectors.
      real(dp), intent(in) :: recLengths(2)

      !> Real-space cutoff radius.
      real(dp), intent(in) :: r_cut

      !> Box counts along the two in-plane directions.
      integer :: nBoxes(2)

      integer :: iDirection

      do iDirection = 1, 2
         nBoxes(iDirection) = floor(1.0_dp/(r_cut*recLengths(iDirection)))
      end do

   end function box_counts_2d

   !> Box a charge falls into along one direction.  The clamp catches a
   !> coordinate that rounds up to exactly one, which would index one box past
   !> the end.
   pure function box_index(frac, nBoxes) result(iBox)

      !> Fractional coordinate along this direction, in [0, 1).
      real(dp), intent(in) :: frac

      !> Number of boxes along this direction.
      integer, intent(in) :: nBoxes

      !> Zero-based box index, in [0, nBoxes-1].
      integer :: iBox

      iBox = min(int(frac*nBoxes), nBoxes - 1)

   end function box_index

   !> Build the decomposition of a three-dimensionally periodic cell.  Charges
   !> are wrapped into the cell as they are binned, so positions may be given
   !> anywhere in space.
   subroutine build_cell_list_3d(positions, nParticle, recVecs, nBoxes, list)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Number of charges to bin.
      integer, intent(in) :: nParticle

      !> Crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(3, 3)

      !> Number of boxes along each direction, from box_counts_3d.
      integer, intent(in) :: nBoxes(3)

      !> The finished decomposition.
      type(TCellList3d), intent(out) :: list

      integer :: iParticle, ix, iy, iz

      list%nBoxes = nBoxes
      allocate (list%frac1(nParticle), list%frac2(nParticle), list%frac3(nParticle))
      allocate (list%next(nParticle))
      allocate (list%head(0:nBoxes(1) - 1, 0:nBoxes(2) - 1, 0:nBoxes(3) - 1))
      list%head = 0
      list%next = 0

      do iParticle = 1, nParticle
         ! Projecting onto the dual basis gives the fractional coordinate;
         ! subtracting its floor wraps it into the cell.
         list%frac1(iParticle) = dot_product(recVecs(1, :), positions(iParticle, :))
         list%frac2(iParticle) = dot_product(recVecs(2, :), positions(iParticle, :))
         list%frac3(iParticle) = dot_product(recVecs(3, :), positions(iParticle, :))
         list%frac1(iParticle) = list%frac1(iParticle) - floor(list%frac1(iParticle))
         list%frac2(iParticle) = list%frac2(iParticle) - floor(list%frac2(iParticle))
         list%frac3(iParticle) = list%frac3(iParticle) - floor(list%frac3(iParticle))

         ix = box_index(list%frac1(iParticle), nBoxes(1))
         iy = box_index(list%frac2(iParticle), nBoxes(2))
         iz = box_index(list%frac3(iParticle), nBoxes(3))

         ! Push the charge onto the front of its box's chain.
         list%next(iParticle) = list%head(ix, iy, iz)
         list%head(ix, iy, iz) = iParticle
      end do

   end subroutine build_cell_list_3d

   !> Build the decomposition of a two-dimensionally periodic cell.  Only the
   !> in-plane coordinates are binned; the coordinate along the open direction
   !> is left to the caller, untouched by the wrapping.
   subroutine build_cell_list_2d(positions, nParticle, recVecs, nBoxes, list)

      !> Cartesian positions, positions(i, :) being the i-th charge.
      real(dp), intent(in) :: positions(:, :)

      !> Number of charges to bin.
      integer, intent(in) :: nParticle

      !> In-plane crystallographic dual basis, one vector per row.
      real(dp), intent(in) :: recVecs(2, 3)

      !> Number of boxes along each in-plane direction, from box_counts_2d.
      integer, intent(in) :: nBoxes(2)

      !> The finished decomposition.
      type(TCellList2d), intent(out) :: list

      integer :: iParticle, ix, iy

      list%nBoxes = nBoxes
      allocate (list%frac1(nParticle), list%frac2(nParticle))
      allocate (list%next(nParticle))
      allocate (list%head(0:nBoxes(1) - 1, 0:nBoxes(2) - 1))
      list%head = 0
      list%next = 0

      do iParticle = 1, nParticle
         list%frac1(iParticle) = dot_product(recVecs(1, :), positions(iParticle, :))
         list%frac2(iParticle) = dot_product(recVecs(2, :), positions(iParticle, :))
         list%frac1(iParticle) = list%frac1(iParticle) - floor(list%frac1(iParticle))
         list%frac2(iParticle) = list%frac2(iParticle) - floor(list%frac2(iParticle))

         ix = box_index(list%frac1(iParticle), nBoxes(1))
         iy = box_index(list%frac2(iParticle), nBoxes(2))

         list%next(iParticle) = list%head(ix, iy)
         list%head(ix, iy) = iParticle
      end do

   end subroutine build_cell_list_2d

end module ewald_cell_list

! This file is part of BandFTN.
! Copyright (C) 2013-2014 Adam Hirst <adam@aphirst.karoo.co.uk>
!
! BandFTN is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! BandFTN is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BandFTN. If not, see <http://www.gnu.org/licenses/>.

module MonkhorstPack
  use Constants
  use Lattice
  use Pseudopotential

  implicit none

  ! container for a Monkhorst-Pack mesh of integer points, pre- and post-symmetrisation
  type Mesh
    type(Latvec), allocatable :: points(:)
    ! scaling factor, for when this will become a list of k-points
    real                      :: factor
    ! the degeneracy (or repeatedness) of each point
    integer,      allocatable :: degen(:)
  contains
    procedure :: Generate => Generate_mesh
    procedure :: Symmetrise => Symmetrise_mesh
  end type Mesh

  ! container for a group of symmetry operations, in matrix representation
  type Symmetry
    ! `iterations` is the number of new points that can be gained from a specific symmetry operation
    integer :: swaps(3), factors(3), iterations
  end type Symmetry

  type(Symmetry), parameter :: symmetries_fcc(14) = &
     [ Symmetry([1,2,3], [ 1, 1, 1], 1), & ! identity, => once unique
       Symmetry([1,3,2], [ 1,-1, 1], 3), & ! rotate about x, y, z by pi/2
       Symmetry([3,2,1], [ 1, 1,-1], 3), & ! 4-fold => thrice unique
       Symmetry([2,1,3], [-1, 1, 1], 3), &
       Symmetry([2,3,1], [ 1, 1, 1], 2), & ! rotate about [111] by 2pi/3
       Symmetry([3,1,2], [ 1,-1,-1], 2), & ! 3-fold => twice unique
       Symmetry([3,1,2], [-1,-1, 1], 2), &
       Symmetry([3,1,2], [-1, 1,-1], 2), &
       Symmetry([2,1,3], [ 1, 1,-1], 1), & ! rotate about [110] by pi
       Symmetry([3,2,1], [ 1,-1, 1], 1), & ! 2-fold => once unique
       Symmetry([1,3,2], [-1, 1, 1], 1), &
       Symmetry([2,1,3], [-1,-1,-1], 1), &
       Symmetry([3,2,1], [-1,-1,-1], 1), &
       Symmetry([1,3,2], [-1,-1,-1], 1) ]

  interface operator (.apply.)
    procedure :: Apply_symmetry
  end interface operator (.apply.)

  private
  public :: Mesh

contains

  pure function Apply_symmetry(left, right) result(transformed)
    ! Returns a mesh-point transformed by a symmetry operation.
    ! The element swapping occurs BEFORE the multiplication factors are applied!
    type(Symmetry), intent(in) :: left
    type(Latvec),   intent(in) :: right
    type(Latvec)               :: transformed

    transformed%hkl = left%factors * right%hkl(left%swaps)
  end function Apply_symmetry

  pure function U(r, q)
    ! Fundamental mesh-generating relation, as in Monkhorst and Pack's paper.
    integer, intent(in) :: r, q
    integer             :: U

    U = (2*r) - q - 1
  end function U

  pure subroutine Generate_mesh(this, q)
    ! Initialises and populates an even mesh of integer points, using the mesh "width" q.
    ! TODO: account for different crystal types (i.e. different RLVs)
    integer,     intent(in)             :: q
    class(Mesh), intent(out)            :: this
    integer,                  parameter :: RLVs(3,3) = reshape([-1,1,1, 1,-1,1, 1,1,-1],[3,3])
    integer                             :: i, j, k

    this%factor = 1.0 / real(2 * q)
    allocate(this%points(q**3))
    allocate(this%degen(q**3))
    do concurrent ( i = 1:q, j = 1:q, k = 1:q )
      ! TODO: pre-compute values of `U(r,q)` for the constant `q`
      ! TODO: populate a 3d array of points (in an efficient order), then reshape into `this%points`
      this%points(i + q*(j-1) + q*q*(k-1)) = Latvec( U(i,q)*RLVs(:,1) + U(j,q)*RLVs(:,2) + U(k,q)*RLVs(:,3) )
    end do
    this%degen(:) = 0
  end subroutine Generate_mesh

  pure subroutine Symmetrise_mesh(this, crystal_type)
    ! Applies crystal symmetry operations to an unsymmetrised mesh, leaving only a symmetry-reduced subset.
    class(Mesh),    intent(in out)              :: this
    integer,        intent(in)                  :: crystal_type
    type(Symmetry),                 allocatable :: symmetries(:)
    integer                                     :: i, j, k, l
    type(Latvec)                                :: current
    integer,                        allocatable :: indices(:)

    ! use the correct set of symmetry operations for the crystal type
    select case (crystal_type)
    case (fcc)
      symmetries = symmetries_fcc
    end select
    ! keep track of the indices (w.r.t. `this%points`) of all the symmetrised points
    allocate(indices(0))
    ! check each meshpoint against all earlier points with degen > 0
    next_point: do i = 1, size(this%points)
      ! using all symmetry operations
      do j = 1, size(symmetries)
        ! make a fresh copy of the current point for each symmetry operation
        current = this%points(i)
        ! for as many times as they uniquely apply
        do k = 1, symmetries(j)%iterations
          ! apply symmetry operation to current point
          current = symmetries(j) .apply. current
          ! check against every previously symmetrised point
          do l = 1, size(indices)
            ! if match, add 1 to the point it matched, and start checking a whole new point
            if ( current == this%points(indices(l)) ) then
              this%degen(indices(l)) = this%degen(indices(l)) + 1
              cycle next_point
            ! also if the inverse matches
            else if ( (-1) * current == this%points(indices(l)) ) then
              this%degen(indices(l)) = this%degen(indices(l)) + 1
              cycle next_point
            end if
          end do
        end do
      end do
      ! if no match, then set its degeneracy to 1, and add to the list of symmetrised points' indices
      this%degen(i) = 1
      indices = [ indices, i ]
    end do next_point
    ! use array-valued index notation to condense the points and their degeneracies into the symmetrised ones only
    ! saves having to perform expensive PACK comparisons and so on
    this%points = this%points(indices)
    this%degen = this%degen(indices)
  end subroutine Symmetrise_mesh

end module MonkhorstPack

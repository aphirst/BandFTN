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
  ! TODO: work out whether this type (and instances of it) can be private to this module
  type Symmetry
    ! `iterations` is the number of new points can be gained from a specific symmetry operation
    integer :: matrix(3,3), iterations
  end type Symmetry

  type(Symmetry), parameter :: fcc_symmetries(14) = &
    (/ Symmetry(reshape([1,0,0,  0,1,0,  0,0,1], [3,3]),1), & ! identity, => once unique
       Symmetry(reshape([1,0,0,  0,0,1,  0,-1,0],[3,3]),3), & ! rotate about x, y, z by pi/2
       Symmetry(reshape([0,0,-1, 0,1,0,  1,0,0], [3,3]),3), & ! 4-fold => thrice unique
       Symmetry(reshape([0,1,0,  -1,0,0, 0,0,1], [3,3]),3), &
       Symmetry(reshape([0,1,0,  0,0,1,  1,0,0], [3,3]),2), & ! rotate about [111] by 2pi/3
       Symmetry(reshape([0,0,-1, 1,0,0,  0,-1,0],[3,3]),2), & ! 3-fold => twice unique
       Symmetry(reshape([0,0,1,  -1,0,0, 0,-1,0],[3,3]),2), &
       Symmetry(reshape([0,0,-1, -1,0,0, 0,1,0], [3,3]),2), &
       Symmetry(reshape([0,1,0,  1,0,0,  0,0,-1],[3,3]),1), & ! rotate about [110] by pi
       Symmetry(reshape([0,0,1,  0,-1,0, 1,0,0], [3,3]),1), & ! 2-fold => once unique
       Symmetry(reshape([-1,0,0, 0,0,1,  0,1,0], [3,3]),1), &
       Symmetry(reshape([0,-1,0, -1,0,0, 0,0,-1],[3,3]),1), &
       Symmetry(reshape([0,0,-1, 0,-1,0, -1,0,0],[3,3]),1), &
       Symmetry(reshape([-1,0,0, 0,0,-1, 0,-1,0],[3,3]),1) /)

  ! container for Density of States data
  ! would benefit from being a parameterised type
  type, extends(Eigencalc) :: DoS
    real,     allocatable :: energies(:), densities(:)
    integer,  allocatable :: populations(:)
    real                  :: dE
  contains
    procedure :: Generate => Generate_DoS
  end type DoS

  interface operator (*)
    procedure :: Apply_symmetry
  end interface operator (*)

  private
  public :: Mesh, Symmetry, fcc_symmetries, DoS, operator(*)

contains

  pure function Apply_symmetry(left, right) result(transformed)
    ! "Applies" a symmetry operation to a point, using matrix multiplication.
    ! overloads the (*) operator
    type(Symmetry), intent(in) :: left
    type(Latvec),   intent(in) :: right
    type(Latvec)               :: transformed

    transformed = Latvec( matmul(left%matrix , right%hkl) )

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
      ! TODO: determine whether it's worth pre-computing the possible outputs of `U`, or working out an alternate vector approach
      this%points(i + q*(j-1) + q*q*(k-1)) = Latvec( U(i,q)*RLVs(:,1) + U(j,q)*RLVs(:,2) + U(k,q)*RLVs(:,3) )
    end do
    this%degen = 0

  end subroutine Generate_mesh

  pure subroutine Symmetrise_mesh(this, symmetries)
    ! Applies crystal symmetry operations to an unsymmetrised mesh, leaving only a symmetry-reduced subset.
    class(Mesh),    intent(in out)              :: this
    type(Symmetry), intent(in)                  :: symmetries(:)
    integer                                     :: i, j, k, l
    type(Latvec)                                :: current
    integer,                        allocatable :: indices(:)

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
          ! apply symmetry operation to current point, note the operator overload
          current = symmetries(j) * current
          ! check against every previously symmetrised point
          do l = 1, size(indices)
            ! if match, add 1 to the point it matched, and start checking a whole new point
            if ( current == this%points(indices(l)) ) then
              this%degen(indices(l)) = this%degen(indices(l)) + 1
              cycle next_point
            ! also if the inverse matches
            else if ( Latvec(-1 * current%hkl) == this%points(indices(l)) ) then
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

  pure subroutine Generate_DoS(this, this_material, magnitude, q, resolution)
    class(DoS),     intent(in out)              :: this
    type(Material), intent(in)                  :: this_material
    integer,        intent(in)                  :: magnitude, q, resolution
    type(Mesh)                                  :: this_mesh
    type(Wavevec),                  allocatable :: kpoints(:)
    integer                                     :: i, j, my_index
    real,                           parameter   :: Emin = -6.0, Emax = 6.0

    ! create integer point-mesh
    call this_mesh%Generate(q)
    call this_mesh%Symmetrise(fcc_symmetries)
    ! convert integer point-mesh into list of k-points
    kpoints = this_mesh%factor * this_mesh%points
    ! compute all eigenenergies at all k-points, this populates this%raw_data
    call this%Compute(kpoints, this_material, magnitude)
    ! allocate arrays
    this%populations = [( 0, i = 1,resolution )]
    allocate(this%densities(resolution))
    ! `resolution` denotes how many population histogram bins to use, within the region [-6.0eV, 6.0eV]
    ! `resolution - 1` is a result of CENTERING the either-end bins at Emin/Emax (rather than placing the bin BOUNDARIES there)
    ! allocate-on-assign
    this%energies = [( Emin + ( (Emax - Emin)*(i-1)/(resolution - 1) ) , i = 1,resolution )]
    ! convert `raw_data` into a population histogram
    ! each set of eigenenergies for a k-point are sorted automatically by LAPACK
    ! so we do each k-point's data (i.e. each column of `raw_data`) separately
    ! TODO: parallelise using an additional array dimension, and a SUM()
    do i = 1, size(this%raw_data,2)
      ! test all values in the column to see whether the population should be incremented
      do j = 1, size(this%raw_data,1)
        my_index = nint( (this%raw_data(j,i) - Emin) * (resolution - 1) / (Emax - Emin) ) + 1
        if (my_index < 1) then
          ! only start doing interesting things once we're in the desired range
          cycle
        else if (my_index > resolution) then
          ! once we find one energy above the desired range, the rest will be too
          exit
        else
          ! add the k-point's degeneracy to the energy subrange's population count
          this%populations(my_index) = this%populations(my_index) + this_mesh%degen(i)
        end if
      end do
    end do
    ! convert the population count into a density using the total (unsymmetrised!) number of kpoints
    this%densities = this%populations / real( q**3 )

  end subroutine Generate_DoS

end module MonkhorstPack

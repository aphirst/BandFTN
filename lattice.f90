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

module Lattice
  use Constants

  ! declaring this at the start of a module means we don't need to do it in every contained routine
  implicit none

  ! container for (reciprocal) lattice vectors
  type Latvec
    integer :: hkl(3)
  contains
    ! ifort complains about multiple PROCEDUREs on the same line (need to check the standard)
    procedure :: Swap
    procedure :: Sort
    procedure :: Shift
    procedure :: Is_valid_RLV
    procedure :: Sym_factor
    procedure :: Asym_factor
  end type Latvec

  ! container for pseudopotential Form-Factors
  type Material
    character(4) :: name
    character(3) :: crystal_type
    ! # of valence electrons; |g''|^2
    integer      :: electrons, gpp_sq(4) = [3,4,8,11]
    ! V(|g''|) symmetric and antisymmetric; lattice constant in Angstroms
    real         :: sym_factors(4), asym_factors(4), A
  contains
    procedure :: Form_factor
  end type Material

  ! defitions of commonly-used materials
  ! TODO: make this work more like a database
  type(Material), parameter :: materials(7) = &
    (/ Material('Ge',   'FCC', 4, [3,4,8,11], [-3.238,0.0,0.0517,0.925], [0.0,0.0,0.0,0.0],       5.65), &
       Material('Si',   'FCC', 4, [3,4,8,11], [-3.048,0.0,0.748,0.98],   [0.0,0.0,0.0,0.0],       5.43), &
       Material('GaAs', 'FCC', 4, [3,4,8,11], [-3.13,0.0,0.136,0.748],   [0.952,0.68,0.0,0.136],  5.64), &
       Material('InP',  'FCC', 4, [3,4,8,11], [-3.606,0.0,0.136,0.816],  [0.952,0.68,0.0,0.136],  5.86), &
       Material('AlAs', 'FCC', 4, [3,4,8,11], [-2.99,0.0,0.585,0.816],   [0.177,0.748,0.0,0.272], 5.66), &
       Material('InAs', 'FCC', 4, [3,4,8,11], [-2.99,0.0,0.0,0.68],      [1.09,0.68,0.0,0.408],   6.04), &
       Material('GaP',  'FCC', 4, [3,4,8,11], [-2.99,0.0,0.408,0.952],   [1.63,0.952,0.0,0.272],  5.44) /)

  ! we want to be able to concisely compare lattice vectors for equivalence
  ! so we overload the == operator
  interface operator (==)
    procedure :: Latvec_equivalence
  end interface operator (==)

  ! we also want to be able to subtract them
  interface operator (-)
    procedure :: Latvec_subtraction
  end interface operator (-)

  private
  public :: Latvec, Material, materials, operator(==), operator(-), Get_groups, Get_signings, Get_permutations

contains

  elemental subroutine Swap(this,i,j)
    ! This really shouldn't need much explanation.
    class(Latvec), intent(in out), target  :: this
    integer,       intent(in)              :: i, j
    integer,                       pointer :: hkl(:)

    ! array-valued index notation
    this%hkl([i,j]) = this%hkl([j,i])

  end subroutine Swap

  elemental subroutine Sort(this)
    ! This is a specialised simple bubble-sorting routine for an {hkl} vector.
    ! The convention is for the elements to be in ascending order.
    ! e.g. {001}, {012}, {223} etc.
    !
    ! Perhaps this could be modified later to match the condensed matter convention
    ! i.e. (1 0 0), (2 -1 0), (-4 3 -3) etc
    ! Sorted in descending element-magnitude, and then by descending sign.
    class(Latvec), intent(in out), target  :: this
    integer,                       pointer :: hkl(:)

    hkl => this%hkl
    do ! until done
      if ( (hkl(1) <= hkl(2)) .and. (hkl(2) <= hkl(3)) ) exit ! done
      if (hkl(1) > hkl(2)) call this%Swap(1,2)
      if (hkl(2) > hkl(3)) call this%Swap(2,3)
    end do

  end subroutine Sort

  elemental subroutine Shift(this, amount)
    ! A type-bound wrapper to make it easier to circle-shift lattice vectors.
    ! Calling this%Shift(n) will shift `this` by n elements, similar to the intrinsic CSHIFT.
    class(Latvec), intent(in out) :: this
    integer,       intent(in)     :: amount

    this%hkl = cshift(this%hkl, amount)

  end subroutine Shift

  elemental function Is_valid_RLV(this)
    ! Determines whether an RLV groups (g) we've generated is "valid" for the crystal type.
    !
    ! g = n1*a + n2*b + n3*c
    ! n = (n1,n2,n3) is a PRIMITIVE lattce vector - components MUST be INTEGERS
    !
    ! Our primitive lattice is FCC, which gives a BCC reciprocal lattice
    ! a=(-1,1,1), b=(1,-1,1), c=(1,1,-1)
    !
    ! Putting this into matrix form and re-arranging gives us a condition that g must satisfy
    ! n = (1/2) M g (note the 1/2)
    ! i.e. for (1/2)matmul(M,g) to be integer, matmul(M,g) must be EVEN.
    !
    ! TODO: perform different checks based on the crystal type
    class(Latvec), intent(in)            :: this
    logical                              :: Is_valid_RLV
    integer,                   parameter :: M(3,3) = reshape( [0,1,1,1,0,1,1,1,0] , [3,3] )
    integer                              :: n(3)

    n = matmul(M, this%hkl)
    Is_valid_RLV = all(modulo(n,2) == 0)

  end function Is_valid_RLV

  elemental function Latvec_equivalence(left, right) result(equiv)
    ! Determines whether two lattice vectors are equivalent.
    ! (overloads the `==` operator)
    type(Latvec), intent(in) :: left, right
    logical                  :: equiv

    equiv = all(left%hkl == right%hkl)

  end function Latvec_equivalence

  elemental function Latvec_subtraction(left, right) result(sub)
    ! Subtracts one lattice vector from another.
    ! (overloads the `-` operator)
    type(Latvec), intent(in) :: left, right
    type(Latvec)             :: sub

    sub = Latvec(left%hkl - right%hkl)

  end function Latvec_subtraction

  pure function Get_groups(M) result(groups)
    ! Generate RLV Miller index 'families' (denoted here as {h k l}) with magnitude sqrt(M),
    ! only keeping the ones that can be validly constructed from integer PRIMITIVE lattice vectors for the crystal type.
    ! This is useful because for the main program we'll want the "expanded" indices IN ORDER OF MAGNITUDE
    ! ("expanded" means that for {0 0 1} we get (0 0 1),(0 1 0),(1 0 0),(0 0 -1),(0 -1 0),(-1 0 0), for instance)
    !
    ! We adopt the convention where h, k and l are in increasing order.
    ! e.g. {0,0,1}, {0,2,2}, {1,2,4} etc.
    !
    ! Optimisation rules for {h,k,l}:
    !
    ! 1) max(max(h),max(k),max(l)) = max(max(l)) = sqrt(M)
    ! 2) min(max(h,k,l)) = min(l) = sqrt(M/n)
    !    where 'n' is the number of yet-undecided elements in {hkl}
    ! 3) Optimisation 2) can be applied on k for each candidate l progressively
    !    with '# of undecided elements' n = 2, AND again on h for each candidate k (n = 1)
    integer,      intent(in)              :: M
    type(Latvec),             allocatable :: groups(:)
    type(Latvec)                          :: group
    integer                               :: h, k, l

    ! for safety's sake we initially allocate `groups` to zero size
    ! also ensures sane function output when no lattice vector groups generate magnitude M
    allocate(groups(0))
    ! count backwards from the upper bound to the lower bound
    ! generates the correct-ordered version of {hkl} first
    ! should make it quicker for later equivalent permutations to be recognised and correctly discarded
    do l = floor(sqrt( real(M) )), floor(sqrt( real(M)/3.0 )), -1
      do k = floor(sqrt( real(M - l**2) )), floor(sqrt( real(M-l**2)/2.0 )), -1
        ! with l and k set, have a condition for what h should be
        h = floor(sqrt( real(M - l**2 - k**2) ))
        ! due to rounding, need to check that this {hkl} gives M
        if ( dot_product([h,k,l],[h,k,l]) == M ) then
          group = Latvec([h,k,l])
          call group%Sort
          ! if not present in the current list of groups, and is a valid RLV, add to the list
          if ( .not. any(groups == group) ) then
            if ( group%Is_valid_RLV() ) groups = [ groups, group ]
          end if
        end if
      end do
    end do

  end function Get_groups

  pure function Get_signings(group) result(signings)
    ! Generate all the unique 'signings' of a (reciprocal) lattice vector group in lexicographic order.
    ! e.g. {1 2 3} gives (-1 -2 -3), (-1 -2 3) etc.
    type(Latvec), intent(in)              :: group
    type(Latvec),             allocatable :: signings(:)
    type(Latvec)                          :: signing
    integer                               :: i, j, k

    allocate(signings(0))
    do i = -1, 1, 2
      do j = -1, 1, 2
        do k = -1, 1, 2
          signing = Latvec([ i*group%hkl(1), j*group%hkl(2), k*group%hkl(3) ])
          call signing%Sort
          if ( .not. any(signings == signing) ) signings = [ signings, signing ]
        end do
      end do
    end do

  end function Get_signings

  pure function Get_permutations(group) result(permutations)
    ! Generate all the unique permutations of an (already 'signed') group in lexicographic order.
    type(Latvec), intent(in)              :: group
    type(Latvec),             allocatable :: permutations(:)
    type(Latvec)                          :: permutation
    integer                               :: i

    allocate(permutations(0))
    permutation = group
    do i = 1, 6
      if (modulo(i,2) == 0) then
        call permutation%Swap(2,3)
      else
        select case (i)
          case (3)
            call permutation%Shift(-1) ! shift right
          case (5)
            call permutation%Shift(1)  ! shift left
        end select
      end if
      if ( .not. any(permutations == permutation) ) permutations = [ permutations, permutation ]
    end do

  end function Get_permutations

  elemental function Sym_factor(this) result(S)
    ! Returns the symmetric structure factor for a value of g'', or g-prime-prime
    ! S(g'') = cos[pi/4 * (n'_x + n'_y + n'_z) ]
    ! where g'' = g - g'
    class(Latvec), intent(in) :: this
    real                      :: S

    S = cos( (pi/4) * real(sum(this%hkl)) )

  end function Sym_factor

  elemental function Asym_factor(this) result(S)
    ! Returns the asymmetric structure factor for a value of g''
    ! S(g'') = sin[pi/4 * (n'_x + n'_y + n'_z) ]
    class(Latvec), intent(in) :: this
    real                      :: S

    S = sin( (pi/4) * real(sum(this%hkl)) )

  end function Asym_factor

  elemental function Form_factor(this, gpp) result(V)
    ! Returns the pseudopotential within a material for a particular g''.
    ! V_(g'') = V_sym(g'') * S_sym(g'') + i V_asym(g'')*S_asym(g'')
    class(Material), intent(in) :: this
    type(Latvec),    intent(in) :: gpp
    complex                     :: V
    integer                     :: i, gpp_sq

    gpp_sq = sum( (gpp%hkl)**2 ) ! |g''|^2
    ! check to see if |g''|^2 matches any of the defined form-factors for the material
    do i = 1, size(this%gpp_sq)
      if (this%gpp_sq(i) == gpp_sq) then
        V = this%sym_factors(i) * gpp%Sym_factor() + (0.0,1.0) * this%asym_factors(i) * gpp%Asym_factor()
        return
      end if
    end do
    ! if it doesn't, return 0
    V = (0.0,0.0)

  end function Form_factor

end module Lattice

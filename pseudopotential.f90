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

module Pseudopotential
  use Constants
  use Lattice

  implicit none

  ! container for wavevectors
  type Wavevec
    real :: k(3)
  end type Wavevec

  ! container for pseudopotential matrix
  type Potmat
    type(Latvec), allocatable :: basis(:)
    integer                   :: N
    ! we use packed storage for the upper triangle of the N-by-N matrix
    complex,      allocatable :: M(:)
    integer,      allocatable :: diag_indices(:)
    ! storing a copy of the lattice constant in here simplifies things somewhat
    real                      :: A
  contains
    procedure :: Create => Create_potmat
    procedure :: EVs => Get_eigenvalues
    ! might be easier to then spread into energy-bands now that eigenvalues are inside here
  end type Potmat

  ! abstract container for an eigen-calculation
  ! allows the bruntwork of "get all eigenvalues for a range of k values" to be separated AND type-bound
  type Eigencalc
    real, allocatable :: raw_data(:,:)
    real              :: Emin, Emax
  contains
    procedure :: Compute => Compute_energies
  end type Eigencalc

  ! we want to be able to multiply a Latvec with integer or real factors to get a Latvec or a Wavevec respectively
  interface operator (*)
    procedure :: Latvec_integer_product
    procedure :: Latvec_real_product
  end interface operator (*)

  ! explicit interface to the packed-storage hermitian LAPACK routine
  interface CHPEV
    subroutine CHPEV(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, INFO)
      character(1), intent(in)     :: JOBZ, UPLO
      integer,      intent(in)     :: N, LDZ
      complex,      intent(in out) :: AP( (N*(N+1))/2 )
      real,         intent(out)    :: W(N)
      complex,      intent(out)    :: Z(LDZ,N), WORK(2*N - 1)
      real,         intent(out)    :: RWORK(3*N - 2)
      integer,      intent(out)    :: INFO
    end subroutine CHPEV
  end interface CHPEV

  private
  public :: Wavevec, Potmat, Eigencalc, operator(*)

contains

  pure subroutine Create_potmat(this, this_material, magnitude)
    ! Initialises the state of type(Potmat) entity, using the properties of `this_material`
    ! and populates its Miller index basis using RLVs of square-magnitude upto 'magnitude'.
    class(Potmat),  intent(in out)              :: this
    type(Material), intent(in)                  :: this_material
    integer,        intent(in)                  :: magnitude
    integer                                     :: i, j
    type(Latvec),                   allocatable :: groups(:), signings(:)
    complex,                        allocatable :: M(:,:)

    ! copy in the lattice constant
    this%A = this_material%A
    ! initialise the RLV basis
    groups = [( Get_groups(i), i=0,magnitude )]
    signings = [( Get_signings(groups(i)), i=1,size(groups) )]
    if (allocated(this%basis)) deallocate(this%basis)
    this%basis = [( Get_permutations(signings(i)), i=1,size(signings) )]
    ! we store the size of the basis inside the type(Potmat) entity for repeated use later
    this%N = size(this%basis)
    ! populate the non-zero non-diagonal elements of a temporary potential matrix
    ! rows are g', columns are g
    ! only need to populate the upper diagonal (since the matrix must have real eigenenergies, it must be Hermitian)
    allocate(M(this%N,this%N))
    do concurrent ( i = 1:size(this%basis), j = 1:size(this%basis), j > i )
      M(i,j) = this_material%Form_factor(this%basis(i) - this%basis(j))
    end do
    ! we cache the indices for the "diagonal" elements of our packed-storage array
    if (allocated(this%diag_indices)) deallocate(this%diag_indices)
    this%diag_indices = [( (i*(i+1))/2 , i=1,this%N )]
    ! we copy the upper diagonal of the temporary matrix into the packed-storage matrix `this%M`
    if (allocated(this%M)) deallocate(this%M)
    this%M = [( M(1:i,i), i = 1,this%N )]
  end subroutine Create_potmat

  function Get_eigenvalues(this, k) result(EVs)
    ! Returns the eigenenergies of a pseudopotential matrix at a given wavevector.
    class(Potmat), intent(in) :: this
    type(Wavevec), intent(in) :: k
    integer                   :: i, INFO
    real                      :: EVs(this%N), RWORK(3*this%N - 2)
    complex                   :: M(size(this%M)), WORK(2*this%N - 1), Z(this%N,this%N)

    ! we make a local copy of the packed pseudopotential matrix
    M(:) = this%M(:)
    ! update all kinetic energy terms, across the matrix's diagonal
    do concurrent ( i = 1:this%N )
      ! technically speaking this should be an addition
      ! for now, the diagonal of the original matrix is always zero, so we can assign without adding
      M(this%diag_indices(i)) = (hc/this%A)**2 * sum( (k%k + real(this%basis(i)%hkl))**2 ) / (2.0 * E_e)
    end do
    ! LAPACK routine for diagonalising complex hermitian matrices in packed storage
    ! 'N' => no eigenvectors
    ! 'U' => the upper diagonal is packed in the input array
    ! eigenvectors would be stored in `Z` if we wanted them
    call CHPEV('N', 'U', this%N, M, EVs, Z, this%N, WORK, RWORK, INFO)
  end function Get_eigenvalues

  subroutine Compute_energies(this, kpoints, this_material, magnitude)
    ! Generates and populates a raw data array of eigenenergies, each column being the results for one k-point.
    class(Eigencalc), intent(in out) :: this
    type(Wavevec),    intent(in)     :: kpoints(:)
    type(Material),   intent(in)     :: this_material
    integer,          intent(in)     :: magnitude
    type(Potmat)                     :: V
    integer                          :: i

    ! initialise the potential / form-factor matrix
    call V%Create(this_material, magnitude)
    allocate(this%raw_data( size(V%basis), size(kpoints) ))
    ! compute eigenenergies for each value of k
    ! this outer loop is PROBABLY parallelisable, depending on the thread-safety of the LAPACK implementation
    ! though if system LAPACK is already multithreaded, there's probably little point setting CONCURRENT here
    do i = 1, size(kpoints)
      this%raw_data(:,i) = V%EVs(kpoints(i))
    end do
    ! renormalise the energy scale based on the number of valence electrons (i.e. number of valence bands)
    ! the top of the last (highest) filled band should be the zero point
    this%raw_data = this%raw_data - maxval( this%raw_data(this_material%electrons,:) )
  end subroutine Compute_energies

  elemental function Latvec_integer_product(left, right) result(prod)
    ! Multiplies an integer scalar by an integer mesh-point (of type Latvec), returning another mesh-point.
    integer,      intent(in) :: left
    type(Latvec), intent(in) :: right
    type(latvec)             :: prod

    prod = Latvec( left * right%hkl )
  end function Latvec_integer_product

  elemental function Latvec_real_product(left, right) result(prod)
    ! Multiplies a real scalar by an integer mesh-point (of type Latvec), returning a real k-point (of type Wavevec).
    real,          intent(in) :: left
    type(Latvec),  intent(in) :: right
    type(Wavevec)             :: prod

    prod = Wavevec( left * right%hkl )
  end function Latvec_real_product

end module Pseudopotential

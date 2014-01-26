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

module DensityofStates
  use Constants
  use Lattice
  use Pseudopotential
  use MonkhorstPack

  implicit none

  ! container for Density of States data
  ! might benefit from being a parameterised type
  type, extends(Eigencalc) :: DoS
    real,     allocatable :: energies(:), densities(:)
  contains
    procedure :: Generate => Generate_DoS
    procedure :: Plot => Plot_DoS
  end type DoS

  private
  public :: DoS

contains

  elemental function Gaussian(x, mu, sigma)
    ! Returns the Gaussian function at `x`, centred at `mu` and with width `sigma`.
    ! Gauss(x) = [ 1 / (sigma *sqrt(2 pi)) ] * exp[ -(x - mu)**2 / (2 * sigma^2) ]
    real, intent(in) :: x, mu, sigma
    real             :: gaussian

    gaussian = ( 1/(sqrt(2*pi) * sigma) ) * exp( -(x - mu)**2 / (2 * sigma**2) )

  end function Gaussian

  subroutine Generate_DoS(this, this_material, magnitude, q, resolution, dE)
    ! Handles the calculation of the density of states for a given material, using a meshsize `q`, and Gaussian smearing.
    class(DoS),     intent(in out)              :: this
    type(Material), intent(in)                  :: this_material
    integer,        intent(in)                  :: magnitude, q, resolution
    real,           intent(in)                  :: dE
    type(Mesh)                                  :: this_mesh
    type(Wavevec),                  allocatable :: kpoints(:)
    integer                                     :: i, j, upper, lower

    ! create integer point-mesh
    call this_mesh%Generate(q)
    call this_mesh%Symmetrise(fcc_symmetries)
    ! convert integer point-mesh into list of k-points
    kpoints = this_mesh%factor * this_mesh%points
    ! compute all eigenenergies at all k-points, this populates this%raw_data
    if (allocated(this%raw_data)) deallocate(this%raw_data)
    call this%Compute(kpoints, this_material, magnitude)
    ! (re-)allocate arrays
    this%energies = [( this%Emin + (this%Emax - this%Emin)*(i-1)/(resolution-1) , i = 1,resolution )]
    this%densities = [( 0.0, i = 1,resolution )]
    ! each set of eigenenergies for a k-point are sorted automatically by LAPACK
    ! so we do each k-point's data (i.e. each column of `raw_data`) separately
    ! TODO: parallelise using an additional array dimension, and a SUM()
    do i = 1, size(this%raw_data,2)
      lower = 1; upper = size(this%raw_data,1)
      ! scan through each column to find the range of indices corresponding to the desired energy range
      do j = 1, size(this%raw_data,1)
        if ( this%raw_data(j,i) < this%Emin ) then
          ! keep updating the `lower` bound until we're inside the range
          ! we use `j+1` instead of `j` since the current `j` is necessarily still outside the range
          lower = j + 1
          cycle
        else if ( this%raw_data(j,i) > this%Emax ) then
          ! once we find one energy above the desired range, the rest will be too
          ! if we're outside the range, the previous value must be the last one inside it (i.e. we use `j-1`)
          upper = j - 1
          exit
        end if
      end do
      ! add the Gaussian-smeared density profiles for this energy range, scaled up by the k-point's degeneracy
      this%densities(:) = this%densities(:) + this_mesh%degen(i) * &
                          [( sum(Gaussian(this%energies(j), this%raw_data(lower:upper,i), dE)) , j = 1,resolution )]
    end do
    ! scale down by the total number of (unsymmetrised!) k-points
    this%densities(:) = this%densities(:) / (q**3)

  end subroutine Generate_DoS

  subroutine Plot_DoS(this)
    class(DoS), intent(in) :: this
  end subroutine Plot_DoS

end module DensityofStates

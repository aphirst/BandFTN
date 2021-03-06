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

module Bandstructure
  use Constants
  use Lattice
  use Pseudopotential

  implicit none

  ! container for a high-symmetry point
  type, extends(Wavevec) :: Highsym
    character(1) :: label
    integer      :: k_index
  end type Highsym

  ! high-symmetry points for the FCC lattice
  type(Highsym), parameter :: sympoints_fcc(6) = [ Highsym([0.5, 0.5, 0.5], 'L', 0), Highsym([0.0, 0.0, 0.0], 'G', 0), &
                                                   Highsym([1.0, 0.0, 0.0], 'X', 0), Highsym([1.0, 0.5, 0.0], 'W', 0), &
                                                   Highsym([0.75,0.75,0.0], 'K', 0), Highsym([0.0, 0.0, 0.0], 'G', 0) ]

  ! container for an energy band
  type Energyband
    real, pointer :: E(:)
  contains
    procedure :: Intersects => Bands_intersect
  end type

  ! container for a bandgap
  type Bandgap
    ! lower- and upper-band numbers are bands(1) and bands(2) respectively
    ! the relevant wavevectors can be found at kpoints(k_indices(1)) and kpoints(k_indices(2))
    integer :: bands(2),  k_indices(2)
    ! E_g is the actual value of the bandgap (in eV), E_offset is the energy at the top of the valence band
    real    :: E_g, E_offset
  end type Bandgap

  ! container for an entire bandstructure, including relevant high-symmetry points
  type, extends(Eigencalc) :: Bandstruc
    type(Highsym),    allocatable :: sympoints(:)
    integer                       :: num_kpoints
    type(Energyband), allocatable :: bands(:)
  contains
    procedure :: Generate => Generate_bandstructure
    procedure :: Plot => Plot_bandstructure
  end type Bandstruc

  private
  public :: Highsym, sympoints_fcc, Bandstruc

contains

  function Expand_highsym(sympoints, resolution) result(kpoints)
    ! Expands a wavevector-array of high-symmetry points into a full basis of wavevectors ("k-points"), for bandstructure plotting.
    ! (linearly interpolates from each high-symmetry point to the next)
    ! `resolution` denotes the number of k-points to generate between two high-symmetry points with unit magnitude difference.
    type(Highsym), intent(in out)              :: sympoints(:)
    integer,       intent(in)                  :: resolution
    type(Wavevec),                 allocatable :: kpoints(:)
    type(Wavevec)                              :: steps(2:size(sympoints))
    integer                                    :: i, j, widths(2:size(sympoints))

    do i = 2, size(sympoints)
      steps(i) = Wavevec(sympoints(i)%k - sympoints(i-1)%k)
      widths(i) = floor( resolution * norm2(steps(i)%k) )
      sympoints(i)%k_index = sympoints(i-1)%k_index + widths(i)
    end do
    allocate(kpoints( sum(widths) ))
    do concurrent ( i = 2:size(sympoints) )
      ! normalise the step-sizes
      steps(i)%k = steps(i)%k / real(widths(i))
      kpoints(sympoints(i-1)%k_index+1:sympoints(i)%k_index) = [( Wavevec(sympoints(i-1)%k+steps(i)%k*real(j)), j=1,widths(i) )]
    end do
  end function Expand_highsym

  pure function Bands_intersect(left, right) result(intersects)
    ! Returns .true./.false. depending on whether the two energybands `left` and `right` intersect.
    class(Energyband), intent(in) :: left, right
    logical                       :: intersects

    ! if left band is below right band
    if ( (minval(left%E) < minval(right%E)) .or. (maxval(left%E) < maxval(right%E)) ) then
      ! if bands intersect
      intersects = ( maxval(left%E) > minval(right%E) )
      ! if right band is below left band
    else
      intersects = ( maxval(right%E) > minval(left%E) )
    end if
  end function Bands_intersect

  subroutine Generate_bandstructure(this, this_material, magnitude, resolution)
    ! Handles the computation of the bandstructure of a given material, using RLVs upto a certain square-magnitude.
    class(Bandstruc), intent(in out),              target :: this
    type(Material),   intent(in)                          :: this_material
    integer,          intent(in)                          :: magnitude, resolution
    type(Wavevec),                    allocatable         :: kpoints(:)
    integer                                               :: i

    select case (this_material%crystal_type)
    case (fcc)
      this%sympoints = sympoints_fcc
    end select
    ! initialise the full basis of wavevectors from the high-symmetry points
    ! also stores the appropriate k-indices of the high-symmetry points into `this%sympoints` as a side-effect
    kpoints = Expand_highsym(this%sympoints, resolution)
    ! we use this to assist later, with plotting
    this%num_kpoints = size(kpoints)
    ! fill `raw_data` with eigenenergies per k-point in each column
    if (allocated(this%raw_data)) deallocate(this%raw_data)
    call this%Compute(kpoints, this_material, magnitude)
    ! allocate as many bands as we care about
    if (allocated(this%bands)) deallocate(this%bands)
    allocate(this%bands( size(this%raw_data,1) ))
    ! point each band to the appropriate row of the raw data
    do concurrent ( i = 1:size(this%bands) )
      this%bands(i)%E => this%raw_data(i,:)
    end do
  end subroutine Generate_bandstructure

  subroutine Plot_bandstructure(this, filename, num_bands)
    ! Outputs the first (most relevant) so-many bands into a datafile, and a gnuplot command file.
    class(Bandstruc), intent(in) :: this
    character(*),     intent(in) :: filename
    integer,          intent(in) :: num_bands
    integer                      :: i

    ! create data file
    open(5, file=trim(filename)//'.bnd', status='replace')
    do i = 1, this%num_kpoints
      ! reading `raw_data` directly by column is vastly more efficient than reading awkward strides from `this%bands`
      write(5,'(i0,*(f10.5))') i, this%raw_data(1:num_bands,i)
    end do
    close(5)

    ! create Gnuplot file
    open(10, file=trim(filename)//'.bnd.gnu', status='replace')
    write(10,'(a)') 'unset key'
    write(10,'(a)') 'set grid'
    write(10,'(a)') 'set ylabel "E [eV]"'
    write(10,'(a,i0,a)') 'set xrange [0:', this%num_kpoints,']'
    ! TODO: dynamic upper-bound for the yrange
    write(10,'(a,f5.1,a,f5.1,a)') 'set yrange [',this%Emin,':',this%Emax,']'
    ! set high-symmetry-points as xtics
    write(10,'(a)',advance='no') 'set xtics ('
    do i = 1,size(this%sympoints)
      write(10,'(3(a),i0)',advance='no') '"', this%sympoints(i)%label, '" ', this%sympoints(i)%k_index
      if (i /= size(this%sympoints)) then
        write(10,'(a)',advance='no') ', '
      else
        write(10,'(a)') ')'
      end if
    end do
    write(10,'(a)') 'set ytics 1'
    write(10,'(a,i0,3(a))') 'plot for [i=2:',num_bands+1,'] "', trim(filename)//'.bnd" u 1:i w l lw 3'
    close(10)
  end subroutine Plot_bandstructure

end module Bandstructure

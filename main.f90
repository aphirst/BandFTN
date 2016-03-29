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

program BandFTN
  use Constants
  use Bandstructure
  use MonkhorstPack

  implicit none

  type(Bandstruc) :: my_bandstructure
  type(Mesh)      :: my_mesh
  integer         :: i

  my_bandstructure%sympoints = sympoints_fcc
  !do i = 1, size(materials)
  !do i = 1,1
    !call my_bandstructure%Generate(materials(i), 50, 10)
    !call my_bandstructure%Plot(materials(i)%name,15)
  !end do
  do i = 1,100
    print *, ''
    print *, i,'=>',i**3
    call my_mesh%Generate(i)
    call my_mesh%Symmetrise(fcc)
  end do

end program BandFTN

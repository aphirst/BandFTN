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

module Constants
  implicit none

  real, parameter :: pi = 4 * atan(1.0)
  real, parameter :: hc = 12398.4193    ! h*c, Planck's constant * speed of light, in (eV)*(Angstroms)
  real, parameter :: E_e = 510998.91013 ! rest mass-energy of an electron, in (eV)
  ! it seems appropriate to use an enumerator for the different possible crystal types
  enum, bind(c)
    enumerator :: fcc, bcc
  end enum

  public

end module Constants

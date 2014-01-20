# This file is part of BandFTN.
# Copyright (C) 2013-2014 Adam Hirst <adam@aphirst.karoo.co.uk>
#
# BandFTN is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BandFTN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BandFTN. If not, see <http://www.gnu.org/licenses/>.

SRCS_f90d1 = \
main.f90 \
constants.f90 \
pseudopotential.f90 \
bandstructure.f90 \
lattice.f90 \
monkhorstpack.f90

OBJS_f90d1 = \
main.o \
constants.o \
pseudopotential.o \
bandstructure.o \
lattice.o \
monkhorstpack.o

SRC_DIR_f90d1 = 
OBJS_DIR = /tmp/BandFTN/obj/
EXE_DIR = /tmp/BandFTN/

EXE = BandFTN
FC = gfortran
IDIR = 
CFLAGS = -O2 -s -march=native -std=f2008 -J$(OBJS_DIR) $(IDIR)
LFLAGS = -llapack 
LIBS = 

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1))

all : $(EXE)

$(EXE) : $(OBJS_f90d1)
	@mkdir -p $(EXE_DIR)
	$(FC) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
main.o: \
    main.f90 \
    bandstructure.o \
    monkhorstpack.o
constants.o: \
    constants.f90
pseudopotential.o: \
    pseudopotential.f90 \
    constants.o \
    lattice.o
bandstructure.o: \
    bandstructure.f90 \
    lattice.o \
    pseudopotential.o
lattice.o: \
    lattice.f90 \
    constants.o
monkhorstpack.o: \
    monkhorstpack.f90 \
    lattice.o


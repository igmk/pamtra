ARCH=$(shell uname -s)



gitHash    := $(shell git show -s --pretty=format:%H)
gitVersion := $(shell git describe)-$(shell git name-rev --name-only HEAD)


FC=gfortran
CC=gcc
FCFLAGS=-c -fPIC -Wuninitialized
ifeq ($(ARCH),Darwin)
	FC=/opt/local/bin/gfortran
	NCFLAGS=-I/opt/local/include/ 
	NCFLAGS_F2PY=-I/opt/local/include/ 
	LFLAGS=-L/opt/local/lib/ -llapack
	LDFLAGS=-lnetcdf -lnetcdff 
else
	NCFLAGS :=  $(shell nc-config --fflags)
	NCFLAGS_F2PY := -I$(shell nc-config --includedir) #f2py does not like -g and -O2
	LFLAGS := -llapack
	LDFLAGS := $(shell nc-config --flibs)
endif




OBJECTS=kinds.o constants.o eps_water.o ref_water.o dielec_water.o


all:  py


%.o: %.f90
	$(FC) $(FCFLAGS) $<

%.o: %.f
	$(FC) $(FCFLAGS) $<


pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	f2py --overwrite-signature -m libRefIndex -h libRefIndex.pyf $(OBJECTS) libRefIndex.f90

py: $(OBJECTS)
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS)  -c --fcompiler=gnu95  libRefIndex.pyf $(OBJECTS) libRefIndex.f90 

pyinstall:
# 	cp *.py ~/lib/python/
	cp *.so ~/lib/python/

ARCH=$(shell uname -s)
F2PY := $(shell which f2py2.7 || echo f2py) #on newer Ubuntu systems, only f2py2.7 is available

src='../../src/'

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

OBJECTS=kinds.o constants.o gasabs_module.o report_module.o


%.o: %.f90
	$(FC) $(FCFLAGS) ${src}$< -o $@

%.f90:
	touch $@

rosen98_gasabs.so: $(OBJECTS)
	$(F2PY) -c --fcompiler=gnu95 rosen98_gasabs.pyf ${src}rosen98_gasabs.f90 $(OBJECTS)

pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	f2py --overwrite-signature -m rosen98_gasabs -h rosen98_gasabs.pyf ${src}rosen98_gasabs.f90 ${src}kinds.f90 ${src}constants.f90 ${src}gasabs_module.f90 ${src}report_module.f90

clean:
	-rm *.so *.mod *.o

install: rosen98_gasabs.so
	cp pamtraGasAbs.py ~/lib/python/
	cp *.so ~/lib/python/


ARCH=$(shell uname -s)
F2PY := $(shell which f2py || echo f2py2.7) #on newer Ubuntu systems, only f2py2.7 is available

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


OBJECTS=kinds.o report_module.o vars_index.o constants.o settings.o viscosity_air.o rho_air.o nan.o rt_utilities.o

%.o: %.f90
	$(FC) $(FCFLAGS) ${src}$< -o $@

%.o: %.f
	$(FC) $(FCFLAGS) ${src}$< -o $@

%.f90:
	touch $@

%.f:
	touch $@

dia2vel.so: $(OBJECTS)
	$(F2PY) $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS)  -c --fcompiler=gnu95  pydia2vel.pyf ${src}dia2vel.f90 $(OBJECTS)
	#$(F2PY) -c --fcompiler=gnu95  pydia2vel.pyf ${src}dia2vel.f90 $(OBJECTS)

pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	$(F2PY) --overwrite-signature -m pyDia2vel -h pydia2vel.pyf ${src}dia2vel.f90 ${src}viscosity_air.f90 ${src}kinds.f90 ${src}report_module.f90 ${src}vars_index.f90 ${src}constants.f90 ${src}settings.f90 ${src}rho_air.f90 ${src}nan.f90 ${src}rt_utilities.f90


install:
	cp *.py ~/lib/python/
	cp *.so ~/lib/python/

clean:
	-rm *.so *.mod *.o

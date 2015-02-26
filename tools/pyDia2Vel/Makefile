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




OBJECTS=kinds.o  report_module.o vars_index.o constants.o settings.o file_mod.o viscosity_air.o rho_air.o 


all:  py


%.o: %.f90
	$(FC) $(FCFLAGS) $<

%.o: %.f
	$(FC) $(FCFLAGS) $<


pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	f2py --overwrite-signature -m pyDia2vel -h pydia2vel.pyf dia2vel.f90 viscosity_air.f90 

py: $(OBJECTS)
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS)  -c --fcompiler=gnu95  pydia2vel.pyf $(OBJECTS) $(pamtraOBJECTS) dia2vel.f90 


install:
	cp *.py ~/lib/python/
	cp *.so ~/lib/python/
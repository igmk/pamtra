ARCH=$(shell uname -s)
F2PY := $(shell which f2py2.7 || echo f2py)

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

OBJECTS=kinds.o settings.o report_module.o eps_water.o constants.o nan.o rt_utilities.o

%.o: %.f90
	$(FC) $(FCFLAGS) ${src}$< -o $@

%.f90:
	touch $@

refWater.so: $(OBJECTS)
	$(F2PY) -c --fcompiler=gnu95  ref_water.pyf ${src}ref_water.f90 $(OBJECTS)

pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	$(F2PY) --overwrite-signature -m refWater -h ref_water.pyf ${src}ref_water.f90 ${src}kinds.f90 ${src}settings.f90 ${src}report_module.f90 ${src}eps_water.f90 ${src}constants.f90 ${src}nan.f90 ${src}rt_utilities.f90

clean:
	-rm *.so *.mod *.o

pyinstall:
# 	cp *.py ~/lib/python/
	cp *.so ~/lib/python/

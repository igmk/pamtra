ARCH=$(shell uname -s)

gitHash    := $(shell git show -s --pretty=format:%H)
gitVersion := $(shell git describe)-$(shell git name-rev --name-only HEAD)

FC=gfortran
CC=gcc
FCFLAGS=-c -fPIC -Wunused -O2 #-I../include/ uninitialized
ifeq ($(ARCH),Darwin)
	FC=/opt/local/bin/gfortran
	NCFLAGS=-I/opt/local/include/ 
	NCFLAGS_F2PY=-I/opt/local/include/ 
	LFLAGS=-L/opt/local/lib/ -llapack -ldfftpack
	LDFLAGS=-lnetcdf -lnetcdff 
else
	NCFLAGS :=  $(shell nc-config --fflags)  -O2
	NCFLAGS_F2PY := -I$(shell nc-config --includedir) #f2py does not like -g and -O2
	LFLAGS := -llapack -ldfftpack
	LDFLAGS := $(shell nc-config --flibs)
endif


#-	NCFLAGS=-I$(HOME)/include
#-	LFLAGS=-L$(HOME)/lib -L/usr/local/lib
#-	LDFLAGS=-lnetcdf -static -lhdf5_hl -lhdf5 -lz -lm -lsz

OBJECTS=kinds.o \
	vars_index.o \
	report_module.o \
	settings.o \
	constants.o \
	nan.o \
	radar_moments.o \
	smooth_savitzky_golay.o \
	radar_hildebrand_sekhon.o \
	dsort.o \
	convolution.o \

BIN=pamtra

all: precompile pamtra py


pamtra: pamtra.f90 $(OBJECTS)
	$(FC) -o $(BIN) pamtra.f90 $(OBJECTS) $(LFLAGS) $(LDFLAGS) 
	mv $(BIN) ../


pamtraDebug: FCFLAGS += -g
pamtraDebug: LFLAGS += -g
pamtraDebug: pamtra
	@echo ""
	@echo "####################################################################################"
	@echo "start debugging with:"
	@echo "gdb ./pamtra"
	@echo "run -n namelist -p profile ..."	
	@echo "####################################################################################"


pamtraProfile: FCFLAGS += -pg
pamtraProfile: LFLAGS += -pg
pamtraProfile: 	pamtra
	@echo ""
	@echo "####################################################################################"
	@echo "analyse gprof output with"
	@echo "gprof ./pamtra | gprof2dot.py | dot -Tpng -o output_old.png"
	@echo "####################################################################################"

pyProfile: NCFLAGS_F2PY += -DF2PY_REPORT_ATEXIT
pyProfile: 	py
	@echo ""
	@echo "####################################################################################"
	@echo "performance report displayed at exit of python"
	@echo "####################################################################################"
pyDebug: NCFLAGS_F2PY += --debug-capi
pyDebug: 	py
	@echo ""
	@echo "####################################################################################"
	@echo "This causes the extension modules to print detailed information while in operation. "
	@echo "####################################################################################"

%.o: %.f90
	$(FC) $(FCFLAGS) $<

%.o: %.f
	$(FC) $(FCFLAGS) $<


pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	f2py --overwrite-signature -m pyRadarMoments -h pyradarmoments.pyf report_module.f90 settings.f90 radar_moments.f90

py: $(OBJECTS)
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS) -c --fcompiler=gnu95  pyradarmoments.pyf $(OBJECTS) 

pyinstall:
	cp *.py ~/lib/python/
	cp *.so ~/lib/python/


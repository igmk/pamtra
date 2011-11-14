ARCH=$(shell uname -s)
FC=gfortran

FCFLAGS=-c -fPIC

gitHash    := $(shell git show -s --pretty=format:%H)
gitVersion := $(shell git describe)-$(shell git name-rev --name-only HEAD)


ifeq ($(ARCH),Darwin)
	FC=/opt/local/bin/gfortran
	NCFLAGS=-I/opt/local/include/ #-g
	LFLAGS=-L/opt/local/lib/
	LDFLAGS=-lnetcdf -lnetcdff
else

 	NCFLAGS=-I/usr/include
 	LFLAGS=-L/usr/lib
	LDFLAGS=-L/usr/lib -lnetcdff -lnetcdf
# 	NCFLAGS=-I$(HOME)/include
# 	LFLAGS=-L$(HOME)/lib -L/usr/local/lib
# 	LDFLAGS=-lnetcdf -static -lhdf5_hl -lhdf5 -lz -lm -lsz


endif




OBJECTS=kinds.o nml_params.o mod_fastem4_coef.o vars_profile.o vars_atmosphere.o mod_io_strings.o constants.o double_moments_module.o\
	vars_output.o conversions.o run_rt3.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o \
	scat_utilities.o dsd_utilities.o mpm93.o get_atmosG0.o get_atmosG.o mie.o mie_densitydep_spheremasseq.o\
	slant_profile.o allocate_profile_vars.o allocate_output_vars.o rosen98_gasabs.o surface.o eps_water.o eps_ice.o eps_mix.o \
	equcom.o land_emis.o equare.o scatdb.o dda_db_liu.o dielec_water.o ref_water.o \
	ref_ice.o e_sat_gg_water.o spec2abs.o interpolation.o deallocate_profile_vars.o deallocate_output_vars.o\
	hydrometeor_extinction.o cloud_ssp.o ice_ssp.o rain_ssp.o snow_ssp.o grau_ssp.o hail_ssp.o legendre2phasefunction.o \
	collect_output.o calculate_active.o fastem4.o specular_surface.o error_msg.o  \
	get_scat_coefs.o dump_profile.o collect_boundary_output.o double_moments.o versionNumber.auto.o

pamtraOBJECTS=formatted_frqstr.o write_nc_results.o

BIN=pamtra

all: precompile pamtra py

	
pamtra: pamtra.f90 $(OBJECTS) $(pamtraOBJECTS)
	$(FC) -o $(BIN) pamtra.f90 $(OBJECTS)  $(pamtraOBJECTS) $(LFLAGS) $(LDFLAGS)
	mv $(BIN) ../
 
scatdb.o:
	gcc -O  -fPIC -c scatdb.c 

%.o: %.f90
	$(FC) $(FCFLAGS) $<
	
write_nc_results.o: write_nc_results.f90
	$(FC) $(FCFLAGS) $(NCFLAGS) write_nc_results.f90
	


pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	f2py -m pyPamtraLib -h pypamtralib.pyf $(pyOBJECTS) pyPamtraLib.f90

py: pyPamtraLib.f90 precompile $(OBJECTS) 
	f2py -c --fcompiler=gnu95 pypamtralib.pyf $(OBJECTS) pyPamtraLib.f90 
	cp pyPamtraLib.so ../py/

pyNode: pyPamtraLib.f90 $(OBJECTS) 
	f2py -c --fcompiler=gnu95 pypamtralib.pyf $(OBJECTS) pyPamtraLib.f90 
	cp pyPamtraLib.so ../py/

pyinstall:
	cp ../py/*.py ~/lib/python/
	cp ../py/*.so ~/lib/python/
	cp ../py/libs/*/*.py ~/lib/python/


precompile: 
	echo "!edit in makefile only!" > versionNumber.auto.f90
	echo "subroutine versionNumber(gitVersion,gitHash)" >> versionNumber.auto.f90
	echo "implicit none" >> versionNumber.auto.f90
	echo "character(40), intent(out) ::gitVersion,gitHash" >> versionNumber.auto.f90
	echo "gitVersion = '$(gitVersion)'" >> versionNumber.auto.f90
	echo "gitHash = '$(gitHash)'" >> versionNumber.auto.f90
	echo "return" >> versionNumber.auto.f90
	echo "end subroutine versionNumber" >> versionNumber.auto.f90
	$(FC) $(FCFLAGS) versionNumber.auto.f90 #otherwise error on first make run!
clean:
	rm -rf $(OBJECTS) $(BIN) *.mod *~

test: test1 test2 test3 test4 cleantest

alltest: test pytest

cleantest:
	@echo "=== Clean up tmp ==="
	rm -fr ../test/tmp/*

test1:
	@echo "=== TEST 1 ==="
	rm -fr ../test/tmp/*
	../pamtra testProfiles.dat ../test/nmls/test1.nml 24 90 150
	../test/compareAscii.sh 1
	@echo "TEST 1 OK"

test2:
	@echo "=== TEST 2 active only === "
	rm -fr ../test/tmp/*
	../pamtra testProfiles.dat ../test/nmls/test2.nml 35
	../test/compareAscii.sh 2
	@echo "TEST 2 OK"

test3:
	@echo "=== TEST 3 passive only === "
	rm -fr ../test/tmp/*
	../pamtra testProfiles.dat ../test/nmls/test3.nml 35
	../test/compareAscii.sh 3
	@echo "TEST 3 OK"

test4:
	@echo "=== TEST 4 Netcdf Data === "
	rm -fr ../test/tmp/*
	../pamtra testProfiles.dat ../test/nmls/test4.nml 35
	../test/compareNetcdf.sh 4
	@echo "TEST 4 OK"

pytest: pyinstall pytest1 pytest2 pytest3 pytest4 cleantest
pytest123: pytest1 pytest2 pytest3 cleantest


pytest1:
	@echo "=== Python TEST 1 ==="
	rm -fr ../test/tmp/*
	python ../test/runPythonTest.py 1
	@echo "TEST 1 OK"

pytest2:
	@echo "=== Python TEST 2 active only === "
	rm -fr ../test/tmp/*
	python ../test/runPythonTest.py 2
	@echo "TEST 2 OK"

pytest3:
	@echo "=== Python TEST 3 passive only === "
	rm -fr ../test/tmp/*
	python ../test/runPythonTest.py 3
	@echo "TEST 3 OK"

pytest4:
	@echo "=== Python TEST 4 Netcdf Data  === "
	rm -fr ../test/tmp/*
	python ../test/runPythonTest.py 4
	../test/compareNetcdf.sh 4
	@echo "TEST 4 OK"


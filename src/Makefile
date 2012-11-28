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


#-	NCFLAGS=-I$(HOME)/include
#-	LFLAGS=-L$(HOME)/lib -L/usr/local/lib
#-	LDFLAGS=-lnetcdf -static -lhdf5_hl -lhdf5 -lz -lm -lsz


OBJECTS=kinds.o constants.o file_mod.o nml_params.o mod_fastem4_coef.o conversions.o double_moments_module.o cosmo_netcdf.o  vars_profile.o \
	vars_atmosphere.o mod_io_strings.o getopt.o parse_options.o radmat.o convolution.o\
	get_atmosphere.o vars_output.o run_rt.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o \
	scat_utilities.o dsd_utilities.o mpm93.o get_atmosG0.o get_atmosG.o mie.o mie_densitydep_spheremasseq.o\
	mie_spec.o radar_spectrum.o radar_simulator.o deallocate_jacobian_vars.o allocate_jacobian_vars.o\
	slant_profile.o allocate_profile_vars.o allocate_output_vars.o rosen98_gasabs.o surface.o eps_water.o eps_ice.o eps_mix.o \
	equcom.o land_emis.o equare.o scatdb.o dda_db_liu.o dielec_water.o ref_water.o ref_ice.o  \
	e_sat_gg_water.o spec2abs.o interpolation.o  rescale_spectra.o deallocate_profile_vars.o deallocate_output_vars.o \
	hydrometeor_extinction_rt3.o cloud_ssp.o ice_ssp.o rain_ssp.o snow_ssp.o grau_ssp.o hail_ssp.o legendre2phasefunction.o \
	collect_output.o calculate_active.o fastem4.o specular_surface.o error_msg.o\
	get_scat_coefs.o dump_profile.o collect_boundary_output.o init_random_seed.o \
	rt4.o radtran4.o radintg4.o radscat4.o tmat_rain_db.o tmat_snow_db.o hydrometeor_extinction_rt4.o scatcnv.o \
	tmatrix_rain.o tmatrix_snow.o matrix_cal.o tmatrix_cal.o lpd.o get_scat_mat.o \
	mass_size_rel.o refractive_index.o avint.o dsort.o rho_air.o viscosity_air.o\
	versionNumber.auto.o dia2vel.o smooth_savitzky_golay.o radar_calc_moments.o \
	radar_hildebrand_sekhon.o

pamtraOBJECTS=formatted_frqstr.o write_nc_results.o 

BIN=pamtra

all: precompile pamtra py


pamtra: pamtra.f90 $(OBJECTS) $(pamtraOBJECTS)
	$(FC) -o $(BIN) pamtra.f90 $(OBJECTS) $(pamtraOBJECTS) $(LFLAGS) $(LDFLAGS) 
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

pyDebug: NCFLAGS_F2PY += --debug-capi
pyDebug: 	py
	@echo ""
	@echo "####################################################################################"
	@echo "This causes the extension modules to print detailed information while in operation. "
	@echo "####################################################################################"


scatdb.o:
	$(CC) -O  -fPIC -c scatdb.c 

%.o: %.f90
	$(FC) $(FCFLAGS) $<

%.o: %.f
	$(FC) $(FCFLAGS) $<

write_nc_results.o: write_nc_results.f90
	$(FC) $(FCFLAGS) $(NCFLAGS) write_nc_results.f90

cosmo_netcdf.o: cosmo_netcdf.f90
	$(FC) $(FCFLAGS) $(NCFLAGS) cosmo_netcdf.f90



pyprecompile: 
	@echo "Make backup before deleting old signature file, auto creating will most likely fail."
	@echo "####################################################################################"
	@echo ""
	f2py -m pyPamtraLib -h pypamtralib.pyf $(pyOBJECTS) pyPamtraLib.f90

py: pyPamtraLib.f90 precompile $(OBJECTS) $(pamtraOBJECTS)
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS)  -c --fcompiler=gnu95  pypamtralib.pyf $(OBJECTS) $(pamtraOBJECTS) pyPamtraLib.f90 
	cp pyPamtraLib.so ../py/

pyNode: pyPamtraLib.f90 $(OBJECTS)  $(pamtraOBJECTS)
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS)  -c --fcompiler=gnu95  pypamtralib.pyf $(OBJECTS) $(pamtraOBJECTS) pyPamtraLib.f90 
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
	rm -rf $(OBJECTS) $(pamtraOBJECTS) $(BIN) *.mod *~

test: test1 test2 test3 test4 cleantest

alltest: test pytest

cleantest:
	@echo "=== Clean up tmp ==="
	rm -fr ../test/tmp/*

test1:
	@echo "=== TEST 1 ==="
	rm -fr ../test/tmp/*
	../pamtra -p testProfiles.dat -n ../test/nmls/test1.nml -f 24,90,150
	../test/compareAscii.sh 1
	@echo "TEST 1 OK"

test2:
	@echo "=== TEST 2 active only === "
	rm -fr ../test/tmp/*
	../pamtra -p testProfiles.dat -n ../test/nmls/test2.nml -f 35
	../test/compareAscii.sh 2
	@echo "TEST 2 OK"

test3:
	@echo "=== TEST 3 passive only === "
	rm -fr ../test/tmp/*
	../pamtra -p testProfiles.dat -n ../test/nmls/test3.nml -f 35
	../test/compareAscii.sh 3
	@echo "TEST 3 OK"

test4:
	@echo "=== TEST 4 Netcdf Data === "
	rm -fr ../test/tmp/*
	../pamtra -p testProfiles.dat -n ../test/nmls/test4.nml -f 35
	../test/compareNetcdf.sh 4
	@echo "TEST 4 OK"

pytest: pyinstall pytest2 pytest3 pytest4 pytest1 cleantest
pytest123: pytest1 pytest2 pytest3 cleantest
pytest234: pytest2 pytest3 pytest4 cleantest

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


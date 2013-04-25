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
	SQLFLAGS := -lfsqlite -lsqlite3 -L../tools/sqlite -L/usr/lib 
else
	NCFLAGS :=  $(shell nc-config --fflags)  -O2
	NCFLAGS_F2PY := -I$(shell nc-config --includedir) #f2py does not like -g and -O2
	LFLAGS := -llapack -ldfftpack
	LDFLAGS := $(shell nc-config --flibs)
	SQLFLAGS := -lfsqlite -lsqlite3 -L../tools/sqlite -L/usr/lib -I../tools/sqlite
endif


#-	NCFLAGS=-I$(HOME)/include
#-	LFLAGS=-L$(HOME)/lib -L/usr/local/lib
#-	LDFLAGS=-lnetcdf -static -lhdf5_hl -lhdf5 -lz -lm -lsz

OBJECTS=kinds.o \
	settings.o \
	report_module.o \
	constants.o \
	mod_fastem4_coef.o \
	gasabs_module.o \
	conversions.o \
	double_moments_module.o \
	cosmo_netcdf.o \
	vars_profile.o \
	vars_atmosphere.o \
	mod_io_strings.o \
	getopt.o \
	parse_options.o \
	rt_utilities.o \
	radmat.o \
	convolution.o \
	get_gasabs.o\
	get_atmosphere.o \
	vars_output.o \
	run_rt.o \
	scat_utilities.o \
	dsd_utilities.o \
	mpm93.o \
	get_atmosG0.o \
	get_surface.o \
	mie.o \
	mie_densitydep_spheremasseq.o\
	dia2vel.o \
	rescale_spectra.o \
	radar_spectrum.o \
	radar_simulator.o \
	deallocate_jacobian_vars.o \
	allocate_jacobian_vars.o\
	slant_profile.o \
	allocate_profile_vars.o \
	allocate_output_vars.o \
	rosen98_gasabs.o \
	surface.o \
	eps_water.o \
	eps_ice.o \
	eps_mix.o \
	equcom.o \
	land_emis.o \
	equare.o \
	scatdb.o \
	dda_db_liu.o \
	dielec_water.o \
	ref_water.o \
	ref_ice.o \
	e_sat_gg_water.o \
	spec2abs.o \
	interpolation.o \
	deallocate_profile_vars.o \
	deallocate_output_vars.o \
	cloud_ssp.o \
	ice_ssp.o \
	rain_ssp.o \
	snow_ssp.o \
	grau_ssp.o \
	hail_ssp.o \
	legendre2phasefunction.o \
	collect_output.o \
	save_active.o \
	fastem4.o \
	specular_surface.o \
	dump_profile.o \
	collect_boundary_output.o \
	random.o \
	rt4.o \
	radtran4.o \
	radintg4.o \
	radscat4.o \
	tmat_rain_db.o \
	tmat_snow_db.o \
	hydrometeor_extinction_rt4.o \
	scatcnv.o \
	tmatrix_rain.o \
	tmatrix_snow.o \
	tmatrix_snow_sql.o \
	tmatrix_ice.o \
	lpq.o \
	get_scat_mat.o \
	mass_size_rel.o \
	refractive_index.o \
	avint.o \
	dsort.o \
	rho_air.o \
	viscosity_air.o\
	versionNumber.auto.o \
	smooth_savitzky_golay.o \
	radar_calc_moments.o \
	radar_hildebrand_sekhon.o \
	write_nc_results.o \
	tmatrix_calc.o \
	amplq.lp.o 
# 	amplq.par.o
#	rt3.o \
	radtran3.o \
	radutil3.o \
	radscat3.o \
	radintg3.o \
	hydrometeor_extinction_rt3.o \
#	tmatrix_cal.o \
#	matrix_cal.o \

BIN=pamtra

all: precompile pamtra py


pamtra: pamtra.f90 $(OBJECTS)
	$(FC) -o $(BIN) pamtra.f90 $(OBJECTS) $(LFLAGS) $(LDFLAGS) $(SQLFLAGS)
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

py: pyPamtraLib.f90 precompile $(OBJECTS)
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS) $(SQLFLAGS) -c --fcompiler=gnu95  pypamtralib.pyf $(OBJECTS) pyPamtraLib.f90 
	cp pyPamtraLib.so ../py/

pyNode: pyPamtraLib.f90 $(OBJECTS) 
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS) $(SQLFLAGS) -c --fcompiler=gnu95  pypamtralib.pyf $(OBJECTS) pyPamtraLib.f90 
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

test: test1 test2 test3 test4 test5 test6 cleantest

alltest: test pytest

cleantest:
	@echo "=== Clean up tmp ==="
	rm -fr ../test/tmp/*

test1:
	@echo "=== TEST 1: no hydrometeors & Rosenkranz ==="
	rm -fr ../test/tmp/*
	../pamtra -p usStandard.dat -n ../test/nmls/test1.nml -f 35,90,150
	../test/compareAscii.sh 1
	@echo "TEST 1 OK"

test2:
	@echo "=== TEST 2: no hydrometeors & Liebe === "
	rm -fr ../test/tmp/*
	../pamtra -p usStandard.dat -n ../test/nmls/test2.nml -f 35,90,150
	../test/compareAscii.sh 2
	@echo "TEST 2 OK"

test3:
	@echo "=== TEST 3 passive only === "
	rm -fr ../test/tmp/*
	../pamtra -p hydrometeors.dat -n ../test/nmls/test3.nml -f 35,90,150
	../test/compareAscii.sh 3
	@echo "TEST 3 OK"

test4:
	@echo "=== TEST 4 active only, simple radar=== "
	rm -fr ../test/tmp/*
	../pamtra -p hydrometeors.dat -n ../test/nmls/test4.nml -f 35,90,150
	../test/compareAscii.sh 4
	@echo "TEST 4 OK"

test5:
	@echo "=== TEST 5 active only, moments === "
	rm -fr ../test/tmp/*
	../pamtra -p hydrometeors.dat -n ../test/nmls/test5.nml -f 35,90,150
	../test/compareAscii.sh 5
	@echo "TEST 5 OK"

test6:
	@echo "=== TEST 6 everything, netcdf === "
	rm -fr ../test/tmp/*
	../pamtra -p hydrometeors.dat -n ../test/nmls/test6.nml -f 35,90,150
	../test/compareNetcdf.sh 6
	@echo "TEST 6 OK"

pytest: pyinstall pytest1 pytest2 pytest3 pytest4 pytest5 pytest6 cleantest
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
	@echo "=== Python TEST 4 active only, simple radar === "
	rm -fr ../test/tmp/*
	python ../test/runPythonTest.py 4
	@echo "TEST 4 OK"
pytest5:
	@echo "=== Python TEST 5 active only, moments === "
	rm -fr ../test/tmp/*
	python ../test/runPythonTest.py 5
	@echo "TEST 5 OK"
pytest6:
	@echo "=== Python TEST 6 everything, netcdf  === "
	rm -fr ../test/tmp/*
	python ../test/runPythonTest.py 6
	../test/compareNetcdf.sh 6
	@echo "TEST 6 OK"


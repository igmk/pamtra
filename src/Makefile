ARCH=$(shell uname -s)
FC=gfortran
FCFLAGS=-c -fdefault-real-8 -fdefault-double-8

gitHash    := $(shell git show -s --pretty=format:%H)
gitVersion := $(shell git describe)


# -g

ifeq ($(ARCH),Darwin)
	NCFLAGS=-I/opt/local/include/ #-g
	LFLAGS=-L/opt/local/lib/
else
	NCFLAGS=-I$(HOME)/include
	LFLAGS=-L$(HOME)/lib
endif

LDFLAGS=-lnetcdf -static

OBJECTS=mod_fastem4_coef.o mod_io_strings.o kinds.o nml_params.o constants.o profile_type.o double_moments_module.o\
	vars_atmosphere.o data.o rt3.o radtran3.o radutil3.o radscat3.o radintg3.o radmat.o \
	scat_utilities.o dsd_utilities.o mpm93.o get_atmosG0.o get_atmosG.o mie.o mie_densitysizedep_spheremasseq.o\
	slant_profile.o allocate_vars_atmosphere.o rosen98_gasabs.o surface.o eps_water.o mie_icefactor.o epsi.o \
	equcom.o land_emis.o equare.o diecon.o emissivity_correction.o ../external/scatdb.o dda_db_liu.o dielec_water.o ref_water.o \
	ref_ice.o e_sat_gg_water.o spec2abs.o\
	hydrometeor_extinction.o cloud_ssp.o ice_ssp.o rain_ssp.o snow_ssp.o grau_ssp.o hail_ssp.o legendre2phasefunction.o \
	vars_output.o collect_output.o calculate_active.o write_nc_results.o fastem4.o specular_surface.o error_msg.o formatted_frqstr.o \
	get_scat_coefs.o dump_profile.o collect_boundary_output.o double_moments.o versionNumber.auto.o

BIN=pamtra

all: precompile pamtra
	
pamtra: pamtra.f90 $(OBJECTS)
	$(FC) -o $(BIN) pamtra.f90 $(OBJECTS) $(LFLAGS) $(LDFLAGS)
	mv $(BIN) ../

precompile: 
	echo "!edit in makefile only!" > versionNumber.auto.f90
	echo "subroutine versionNumber(gitVersion,gitHash)" >> versionNumber.auto.f90
	echo "implicit none" >> versionNumber.auto.f90
	echo "character(40), intent(out) ::gitVersion,gitHash" >> versionNumber.auto.f90
	echo "gitVersion = '$(gitVersion)'" >> versionNumber.auto.f90
	echo "gitHash = '$(gitHash)'" >> versionNumber.auto.f90
	echo "return" >> versionNumber.auto.f90
	echo "end subroutine versionNumber" >> versionNumber.auto.f90


%.o: %.f90
	$(FC) $(FCFLAGS) $<

scatdb.o: scatdb.c
	gcc -c ../external/scatdb.c -o ../external/scatdb.o
	
dda_db_liu.o: dda_db_liu.f90
	$(FC) $(FCFLAGS) dda_db_liu.f90 ../external/scatdb.o

write_nc_results.o: write_nc_results.f90
	$(FC) $(FCFLAGS) $(NCFLAGS) write_nc_results.f90
	
clean:
	rm -rf $(OBJECTS) $(BIN) *.mod *~

test:
	@echo "=== TEST 1 ==="
	../pamtra rt_comp_trio.dat ../test/nmls/test1.nml 24 90 150
	cat ../test/tmp/* > ../test/testResult.txt
	diff ../test/testResult.txt ../test/referenceOutput/test1.txt
	rm ../test/testResult.txt ../test/tmp/*

	@echo "=== TEST 2 active only === "
	../pamtra rt_comp_trio.dat ../test/nmls/test2.nml 35
	cat ../test/tmp/* > ../test/testResult.txt
	diff ../test/testResult.txt ../test/referenceOutput/test2.txt
	rm ../test/testResult.txt ../test/tmp/*

	@echo "=== TEST 3 passive only === "
	../pamtra rt_comp_trio.dat ../test/nmls/test3.nml 35
	cat ../test/tmp/* > ../test/testResult.txt
	diff ../test/testResult.txt ../test/referenceOutput/test3.txt
	rm ../test/testResult.txt ../test/tmp/*



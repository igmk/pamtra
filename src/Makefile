ARCH=$(shell uname -s)
HOME=/home/mech
FC=gfortran
FCFLAGS=-c -fdefault-real-8 -fdefault-double-8

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
	get_scat_coefs.o dump_profile.o collect_boundary_output.o double_moments.o

BIN=pamtra

all: pamtra.f90 $(OBJECTS)
	$(FC) -o $(BIN) pamtra.f90 $(OBJECTS) $(LFLAGS) $(LDFLAGS)
	mv $(BIN) ../
	
pamtra: pamtra.f90 $(OBJECTS)
	$(FC) -o $(BIN) pamtra.f90 $(OBJECTS) $(LFLAGS) $(LDFLAGS)
	mv $(BIN) ../

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

ARCH=$(shell uname -s)

gitHash    := $(shell git show -s --pretty=format:%H)
gitVersion := $(shell git describe)-$(shell git name-rev --name-only HEAD)

FC=gfortran
CC=gcc
FCFLAGS=-c -fPIC -Wunused -O2 -cpp #-I../include/ uninitialized
ifeq ($(ARCH),Darwin)
	FC=/opt/local/bin/gfortran-mp-4.8
	NCFLAGS=-I/opt/local/include/ 
	NCFLAGS_F2PY=-I/opt/local/include/ 
	LFLAGS= -L../lib -ldfftpack  -L/opt/local/lib/ -llapack
	LDFLAGS=-lnetcdf -lnetcdff  -lz
else
	NCFLAGS :=  $(shell nc-config --fflags)  -O2
	NCFLAGS_F2PY := -I$(shell nc-config --includedir) #f2py does not like -g and -O2
	LFLAGS := -llapack -L../lib -ldfftpack
	LDFLAGS := $(shell nc-config --flibs) -lz
endif




OBJECTS=kinds.o \
        vars_index.o \
	report_module.o \
	settings.o \
	constants.o \
	nan.o \
	zlib_stuff.o \
	mod_fastem4_coef.o \
	gasabs_module.o \
	conversions.o \
	descriptor_file.o \
	vars_atmosphere.o \
	vars_rt.o \
        vars_hydroFullSpec.o \
	mod_io_strings.o \
	getopt.o \
	parse_options.o \
	rt_utilities.o \
	radmat.o \
	convolution.o \
	get_gasabs.o\
	vars_output.o \
	run_rt.o \
	scat_utilities.o \
	mpm93.o \
	eps_water.o \
	get_surface.o \
	mie_scat_utilities.o \
	mie_spheres.o \
	dia2vel.o \
	rescale_spectra.o \
	radar_moments.o \
	radar_spectrum.o \
	radar_simulator.o \
	rosen98_gasabs.o \
	surface.o \
	eps_ice.o \
	eps_mix.o \
	equcom.o \
	land_emis.o \
	equare.o \
	scatdb.o \
	ref_water.o \
	ref_ice.o \
	e_sat_gg_water.o \
	interpolation.o \
	collect_output.o \
	save_active.o \
	fastem4.o \
	specular_surface.o \
	random.o \
	rt4.o \
	radtran4.o \
	radintg4.o \
	radscat4.o \
	drop_size_dist.o \
	make_dist.o \
	make_dist_param.o \
	make_mass_size.o \
	calc_moment.o \
	make_soft_spheroid.o \
	check_print.o \
	tmatrix.o \
	scatProperties.o \
	hydrometeor_extinction.o \
	scatcnv.o \
	tmatrix_lpq.o \
	get_scat_mat.o \
	refractive_index.o \
	dsort.o \
	rho_air.o \
	viscosity_air.o\
	versionNumber.auto.o \
	smooth_savitzky_golay.o \
	radar_hildebrand_sekhon.o \
	write_nc_results.o \
	tmatrix_amplq.lp.o \
	deallocate_everything.o

BIN=pamtra

all: dfftpack pamtra py py_usStandard

dfftpack:
	cd ../tools/dfftpack && $(MAKE) && cp libdfftpack.a ../../lib/libdfftpack.a

pamtra: pamtra.f90 precompile $(OBJECTS)
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
	@echo "Make backup before deleting old signature file, auto creating can fail."
	@echo "#######################################################################"
	@echo ""
	f2py --overwrite-signature -m pyPamtraLib -h pypamtralib.pyf report_module.f90 deallocate_everything.f90 vars_output.f90 vars_atmosphere.f90 settings.f90 descriptor_file.f90 vars_hydroFullSpec.f90 pyPamtraLib.f90

py: pyPamtraLib.f90 precompile $(OBJECTS)
	f2py $(NCFLAGS_F2PY) $(LDFLAGS) $(LFLAGS) -c --fcompiler=gnu95  pypamtralib.pyf $(OBJECTS) pyPamtraLib.f90 
	cp pyPamtraLib.so ../pyPamtra/

py_usStandard:
	cd ../tools/py_usStandard/ && $(MAKE) all


pyinstall:
	cp -r ../pyPamtra ~/lib/python/
	cd ../tools/py_usStandard/ && $(MAKE) install


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
	cd ../tools/dfftpack/ && $(MAKE) clean
	cd ../tools/py_usStandard/ && $(MAKE) clean


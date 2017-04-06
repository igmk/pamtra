ARCH=$(shell uname -s)

#for some strange reasons, python modules compiled with f2py do not work when using a separate build directory
OBJDIR := src/
SRCDIR := src/
BINDIR := bin/
LIBDIR := lib/
PYTDIR := python/pyPamtra
PYINSTDIR := ~/lib/python/

gitHash    := $(shell git show -s --pretty=format:%H)
gitVersion := $(shell git describe)-$(shell git name-rev --name-only HEAD)

NCCONF = $(shell which nf-config || which nc-config) # on newer Ubuntu version C and Fortran libraries have their own configure scripts
F2PY := $(shell which f2py2.7 || which f2py) # on newer Ubuntu systems, only f2py2.7 is available
FC=gfortran
CC=gcc
FCFLAGS=-c -fPIC -Wunused  -cpp -J$(OBJDIR) -I$(OBJDIR)
#FCFLAGS=-g -c -fPIC -Wunused -O0 -cpp -J$(OBJDIR) -I$(OBJDIR)

NCFLAGS :=  $(shell $(NCCONF) --fflags) 
LFLAGS := -L/usr/lib/ -llapack -L$(LIBDIR) -L../$(LIBDIR) -ldfftpack -lblas -lz
LDFLAGS := $(shell $(NCCONF) --flibs) 
# it's messi but needed for ubuntu 16.04
to_remove:=-Wl,-Bsymbolic-functions -Wl,-z,relro
LDFLAGS := $(subst $(to_remove),,$(LDFLAGS))




OBJECTS=kinds.o \
	vars_index.o \
	report_module.o \
	rt_utilities.o \
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
	radmat.o \
	convolution.o \
	get_gasabs.o\
	vars_output.o \
	azimuth_emissivity_module.o \
	hyperbolic_step.o \
	slope_variance.o \
	reflection_correction_module.o \
	large_scale_correction_module.o \
	small_scale_correction_module.o \
	foam_utility_module.o \
	liu.o \
	fresnel.o \
	fastemx.o \
	ocean_sfc_optics.o \
	land_sfc_optics.o \
	sfc_optics.o \
	sfc_matrices.o \
	run_rt.o \
	scat_utilities.o \
	mpm93.o \
	eps_water.o \
	mie_scat_utilities.o \
	mie_spheres.o \
	scatdb.o \
	dda_db_liu.o \
	dda_db_hong.o \
	hongdb.o \
	dia2vel.o \
	rescale_spectra.o \
	radar_moments.o \
	radar_spectrum.o \
	radar_spectral_broadening.o \
	radar_simulator.o \
	rosen98_gasabs.o \
	surface.o \
	eps_ice.o \
	eps_mix.o \
	equcom.o \
	land_emis.o \
	equare.o \
	ref_water.o \
	ref_ice.o \
	e_sat_gg_water.o \
	interpolation.o \
	collect_output.o \
	save_active.o \
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
	rayleigh_gans.o \
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
	tmatrix_amplq.lp.o \
	deallocate_everything.o
OBJECTS_NC=write_nc_results.o


FOBJECTS=$(addprefix $(OBJDIR),$(OBJECTS))
FOBJECTS_NC=$(addprefix $(OBJDIR),$(OBJECTS_NC))

BIN=pamtra

all: dfftpack pamtra py py_usStandard 

warning:
ifndef PAMTRA_DATADIR
	@echo "########################################################"
	@echo "Do not forget to define PAMTRA_DATADIR environment variable"
	@echo "########################################################"
endif

print-%  : ; @echo $* = $($*)

dfftpack: | $(LIBDIR)
	cd tools/dfftpack && $(MAKE)
	cp tools/dfftpack/libdfftpack.a $(LIBDIR)

pamtra: FCFLAGS += -O2
pamtra: NCFLAGS += -O2
pamtra: warning dfftpack $(FOBJECTS) $(BINDIR)$(BIN) | $(BINDIR) 

$(OBJDIR)versionNumber.auto.o: .git/HEAD .git/index
	echo "!edit in makefile only!" > $(SRCDIR)versionNumber.auto.f90
	echo "subroutine versionNumber(gitVersion,gitHash)" >> $(SRCDIR)versionNumber.auto.f90
	echo "implicit none" >> $(SRCDIR)versionNumber.auto.f90
	echo "character(40), intent(out) ::gitVersion,gitHash" >> $(SRCDIR)versionNumber.auto.f90
	echo "gitVersion = '$(gitVersion)'" >> $(SRCDIR)versionNumber.auto.f90
	echo "gitHash = '$(gitHash)'" >> $(SRCDIR)versionNumber.auto.f90
	echo "return" >> $(SRCDIR)versionNumber.auto.f90
	echo "end subroutine versionNumber" >> $(SRCDIR)versionNumber.auto.f90
	$(FC) $(FCFLAGS) $(SRCDIR)versionNumber.auto.f90 -o $(OBJDIR)versionNumber.auto.o #otherwise error on first make run!

$(OBJDIR):
	mkdir -p $(OBJDIR)
$(LIBDIR):
	mkdir -p $(LIBDIR)
$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)$(BIN): $(FOBJECTS) $(FOBJECTS_NC) | $(BINDIR)
	$(FC) -I$(OBJDIR) -o $(BINDIR)$(BIN) $(SRCDIR)pamtra.f90 $(FOBJECTS) $(FOBJECTS_NC) $(LFLAGS) $(LDFLAGS)

$(OBJDIR)scatdb.o:  $(SRCDIR)scatdb.c  | $(OBJDIR)
	$(CC) -O  -fPIC -c $< -o $@

$(OBJDIR)%.o:  $(SRCDIR)%.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS) $< -o $@

$(OBJDIR)%.o:  $(SRCDIR)%.f | $(OBJDIR)
	$(FC) $(FCFLAGS) $< -o $@

$(OBJDIR)write_nc_results.o:  $(SRCDIR)write_nc_results.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS) $(NCFLAGS) $< -o $@



pamtraDebug: FCFLAGS += -g -fbacktrace -fbounds-check
pamtraDebug: LFLAGS += -g -fbacktrace -fbounds-check
pamtraDebug: pamtra
	@echo ""
	@echo "####################################################################################"
	@echo "start debugging with:"
	@echo "gdb ./pamtra"
	@echo "run -n namelist -p profile ..."
	@echo "or with valgrind:"
	@echo "valgrind --leak-check=yes ./pamtra ..."
	@echo "####################################################################################"


pamtraProfile: FCFLAGS += -pg
pamtraProfile: LFLAGS += -pg
pamtraProfile: 	pamtra
	@echo ""
	@echo "####################################################################################"
	@echo "analyse gprof output with"
	@echo "gprof ./pamtra | gprof2dot.py | dot -Tpng -o output_old.png"
	@echo "####################################################################################"

pyProfile: LFLAGS += -DF2PY_REPORT_ATEXIT
pyProfile: 	py
	@echo ""
	@echo "####################################################################################"
	@echo "performance report displayed at exit of python"
	@echo "####################################################################################"

pyDebug: LFLAGS += -g -fbacktrace -fbounds-check
pyDebug: LFLAGS += --debug-capi
pyDebug: 	py
	@echo ""
	@echo "####################################################################################"
	@echo "This causes the extension modules to print detailed information while in operation. "
	@echo "####################################################################################"



$(OBJDIR)pypamtralib.pyf:  $(FOBJECTS)
	@echo "####################################################################################"
	@echo "Note there is a bug in numpy 1.10.1, intent in or out is not recognized"
	@echo "Note there is a bug in numpy 1.12.0, length of arrays is not recognized by f2py"
	@echo "####################################################################################"
	$(F2PY) --overwrite-signature -m pyPamtraLib -h $(OBJDIR)pypamtralib.pyf $(SRCDIR)report_module.f90 $(SRCDIR)vars_index.f90 $(SRCDIR)viscosity_air.f90 $(SRCDIR)convolution.f90 $(SRCDIR)deallocate_everything.f90 $(SRCDIR)vars_output.f90 $(SRCDIR)vars_atmosphere.f90 $(SRCDIR)settings.f90 $(SRCDIR)descriptor_file.f90 $(SRCDIR)vars_hydroFullSpec.f90 $(SRCDIR)radar_moments.f90 $(SRCDIR)eps_water.f90  $(SRCDIR)radar_hildebrand_sekhon.f90 $(SRCDIR)dia2vel.f90 $(SRCDIR)pyPamtraLib.f90

py: FCFLAGS += -O2
py: NCFLAGS += -O2
py: $(PYTDIR)pyPamtraLib.so

$(PYTDIR)pyPamtraLib.so:  $(SRCDIR)pyPamtraLib.f90 $(OBJDIR)pypamtralib.pyf $(FOBJECTS) | $(BINDIR)
	cd $(OBJDIR) && $(F2PY) $(LFLAGS) -c --fcompiler=gnu95  ../$(OBJDIR)pypamtralib.pyf $(OBJECTS) ../$(SRCDIR)pyPamtraLib.f90
	mv $(OBJDIR)/pyPamtraLib.so $(PYTDIR)
	cp $(PYTDIR)/pamtra.py $(BINDIR)


py_usStandard:
	cd tools/py_usStandard/ && $(MAKE) all

pyinstall: warning dfftpack py py_usStandard
	mkdir -p $(PYINSTDIR)
	cp -r $(PYTDIR) $(PYINSTDIR)
	cd tools/py_usStandard/ && $(MAKE) install
	cp tools/pyRadarMoments/radarMoments.py	 $(PYINSTDIR)

clean:
	-rm -f $(OBJDIR)*.o
	-rm -f $(OBJDIR)*.mod
	-rm -f $(BINDIR)pamtra*
	cd tools/dfftpack/ && $(MAKE) clean
	cd tools/py_usStandard/ && $(MAKE) clean

htmldoc:
	cd doc && $(MAKE) html


all: usStandardAtmosphere.so

usStandardAtmosphere.so: usStandard.f90 
	f2py2.7 -c --fcompiler=gnu95  usStandardAtmosphere.pyf usStandard.f90 


install: usStandardAtmosphere.so
	cp *.py ~/lib/python/
	cp *.so ~/lib/python/

clean: 
	-rm *.so

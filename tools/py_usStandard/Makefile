
all: compile

compile: 
	f2py -c --fcompiler=gnu95  usStandardAtmosphere.pyf usStandard.f90 


install:
	cp *.py ~/lib/python/
	cp *.so ~/lib/python/

clean: 
	rm *.so

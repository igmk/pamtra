F2PY := $(shell which f2py2.7 || echo f2py) #on newer Ubuntu systems, only f2py2.7 is available

all: usStandardAtmosphere.so

usStandardAtmosphere.so: usStandard.f90
	$(F2PY) -c --fcompiler=gnu95  usStandardAtmosphere.pyf usStandard.f90


install: usStandardAtmosphere.so
	cp *.py ~/lib/python/
	cp *.so ~/lib/python/

clean:
	-rm *.so

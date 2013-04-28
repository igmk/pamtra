# lib and makefile taken from flibs: http://flibs.sourceforge.net/fsqlite.html
#
# Makefile for "fsqlite" library
#
# NOTE:
# You need to adjust the macro DEFINES in this makefile:
# - To select the calling conventions for Fortran and C
# - To set the directory where the Sqlite3 header file can be found
#
# The calling conventions for the Fortran compiler need to
# match those for the C compiler. To this end the following
# macros have been defined in the source code:
#
# - LOWERCASE:
#   Define this macro if the Fortran compiler translates the names to
#   lower case and adds a single underscore
#   E.g.
#   A subroutine SQLITE3_OPEN in the Fortran code becomes sqlite3_open_
#   in the object file.
#
# - DBL_UNDERSCORE:
#   Define this macro if the Fortran compiler translates the names to
#   lower case but adds two underscores if the name itself contains
#   an underscore.
#   E.g.
#   A subroutine SQLITE3_OPEN in the Fortran code becomes sqlite3_open__
#   in the object file.
#
# - WIN32:
#   Under Windows the decoration is often more complicated. Most
#   compilers under Windows define this macro automatically, so it is
#   used to identify this type of decoration, but there can be
#   confusion: under Cygwin and MingW this macro is also defined, but
#   the decorations are not necessarily those of "__stdcall".
#   (Using either of the two others overrides this).
#
# Note:
# Things are further complicated, under Windows at least, if you
# want to create a DLL. The C and Fortran source code as given here
# is not enough then (you need to export routines explicitly from the
# DLL).
#
# $Id: makefile,v 1.1 2008/05/04 13:23:56 arjenmarkus Exp $
#
SPACE	=	\

SEP	=	/

FC	=	gfortran
FFLAGS_NORMAL	=-fPIC -Wunused -O3 
FFLAGS_DEBUG	=	-c -g
FFLAGS_OPTIMISE	=	-c -O

LD	=	gfortran
LDFLAGS_NORMAL	=
LDFLAGS_DEBUG	=	-g
LDFLAGS_OPTIMISE	=	
#LDOUTPUT	=	-o$(SPACE)
LDOUTPUT	=	-o $@

LIB	=	ar r

OBJEXT	=	.o
EXEEXT	=	
MODEXT	=	.mod

DELETE	=	rm -f
FFLAGS	=	$(FFLAGS_NORMAL)
LDFLAGS	=	$(LDFLAGS_NORMAL)
CFLAGS =	-fPIC -Wunused -O3
LIBRARIES	=	libfsqlite.a \

#
# Use the correct macro for interfacing between Fortran and C and
# add the directory where the Sqlite3 header file can be found
#
#DEFINES	=	-DDBL_UNDERSCORE -I/usr/include
DEFINES	=	-DLOWERCASE  -I/usr/include
# DEFINES	=	-DLOWERCASE  -I/home/mmaahn/Downloads/sqlite-autoconf-3071602/

#DEFINES	=	-I../../../sqlite


all:	$(LIBRARIES)

libfsqlite.a:	fsqlite.o csqlite.o
	ar r libfsqlite.a fsqlite.o csqlite.o

fsqlite.o:	fsqlite.f90
	$(FC) -c $(FFLAGS) fsqlite.f90

csqlite.o:	csqlite.c
	$(CC) -c $(CFLAGS) csqlite.c $(DEFINES)

clean:	
	$(DELETE) *$(OBJEXT)
	$(DELETE) $(PROGRAMS)
	$(DELETE) *$(MODEXT)
	$(DELETE) createTestFile runTest somedata.* *.a

example:
	gfortran -fPIC -Wunused -O1 -o runTest  csvtable.f90 -lfsqlite -lsqlite3 -L. -L/usr/lib   
	gfortran -fPIC -Wunused -O1  csvdata.f90  -o createTestFile
# 

test:
	./createTestFile
	./runTest
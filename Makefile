BASEPATH=/Users/hzl/plasma/PlasmaDischarged/myPIC/Final_Test
CXX=clang++
CFLAGS=-std=c++17 -Wall -g -O2
PROG=main

OBJS=main.o espic_math.o espic_info.o parse.o str_split.o \
     mesh.o param_particle.o species.o particles.o ambient.o \
     tile.o reaction.o cross_section.o collision.o
	
EIGEN_PATH=${BASEPATH}/ThirdParty
EIGEN=${EIGEN_PATH}/Eigen3.3.7

INCLUDES=-I$(EIGEN)

LIBOBJDIR=Object 
LIBOBJ=libobject.a
LIBINJDIR=Inject
LIBINJ=libinject.a
LIBS=-L$(LIBINJDIR) -L$(LIBOBJDIR)
LINKOPTS=-linject -lobject -lm

all : libinject libobject $(OBJS)
	$(CXX) $(CFLAGS) $(LIBS) $(OBJS) -o $(PROG) $(LINKOPTS)

libinject :
	@cd $(LIBINJDIR); \
	make -f Makefile.inject $(LIBINJ) \
	CXX='$(CXX)' CFLAGS='$(CFLAGS)' \
	INCLUDES='$(INCLUDES)' LIBINJ='$(LIBINJ)'

libobject :
	@cd $(LIBOBJDIR); \
	make -f Makefile.object $(LIBOBJ) \
	CXX='$(CXX)' CFLAGS='$(CFLAGS)' \
	INCLUDES='$(INCLUDES)' LIBINJ='$(LIBOBJ)'

.cpp.o :
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

run :
	./$(PROG)

inject_clean :
	cd $(LIBINJDIR); make -f Makefile.inject clean

object_clean :
	cd $(LIBOBJDIR); make -f Makefile.object clean

clean : inject_clean object_clean
	/bin/rm -f *.o

distclean: clean

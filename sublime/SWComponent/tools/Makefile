CC  = sw5cc
CXX = mpiCC
FC  = mpiswgfortran
AR  = swar
RANLIB = swranlib

ARFLAGS  = cru
CFLAGS   = -g -O2 -OPT:IEEE_arithmetic=2
CHOSTFLAGS   = ${CFLAGS} -host
CSLAVEFLAGS   = ${CFLAGS} -slave
CXXFLAGS = -g -O2 -OPT:IEEE_arithmetic=2
FFLAGS = -g -O2
LDFLAGS  = 

CPPFiles = timers/Timers.cxx
FFiles = timers/Timers_wrapper.f90

INCLUDE= -I../include/ \

CPPOBJS=${CPPFiles:.cxx=.o}
FOBJS= ${FFiles:.f90=.o}

libtools='../lib/libtools.a'

default:${libtools}

${libtools}:${CPPOBJS} ${FOBJS}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@
${OBJS}: %.o:%.cpp
	${CXX} ${CXXFALGS} ${INCLUDE} -c $< -o $@
${HOSTOBJS}: %.o:%.c 
	${CC} ${CHOSTFLAGS} ${INCLUDE} -c $< -o $@
${SLAVEOBJS}: %.o:%.c
	${CC} ${CSLAVEFLAGS} ${INCLUDE} -c $< -o $@
${CPPOBJS}: %.o:%.cxx
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $< -o $@
${FOBJS}: %.o:%.f90
	${FC} ${FFLAGS} ${INCLUDE} -c $< -o $@

.PHONY:clean
clean:
	rm -f ${CPPOBS} ${FOBJS} ${libtools}

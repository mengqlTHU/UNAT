CC  = sw5cc
CXX = mpiCC
AR  = swar
RANLIB = swranlib

ARFLAGS  = cru
CFLAGS   = -O2 -OPT:IEEE_arithmetic=2
CHOSTFLAGS   = ${CFLAGS} -host
CSLAVEFLAGS   = ${CFLAGS} -slave
CXXFALGS = -O2 -OPT:IEEE_arithmetic=2
LDFLAGS  = 

files = gra_gg_pre/gra_gg_pre_host.cpp upw_conv/upw_conv_host.cpp \
	rhs_vis/rhs_vis_host.cpp gra_gg/gra_gg_host.cpp uc2u/uc2u_host.cpp \
	lhs_scalar/lhs_scalar_host.cpp exp_rk43/exp_rk43_host.cpp \
	rhs_rk3/rhs_rk3_host.cpp imp_LHS/imp_LHS_host.cpp
CFiles =
ptrFiles = upw_conv/upw_conv_ptr.c rhs_vis/rhs_vis_ptr.c \
	gra_gg_pre/gra_gg_pre_ptr.c gra_gg/gra_gg_ptr.c uc2u/uc2u_ptr.c \
	lhs_scalar/lhs_scalar_ptr.c exp_rk43/exp_rk43_ptr.c \
	rhs_rk3/rhs_rk3_ptr.c imp_LHS/imp_LHS_ptr.c

INCLUDE= -I../include/ \

OBJS=${files:.cpp=.o}
HOSTOBJS=${CFiles:.c=.o}
SLAVEOBJS=${CSlaveFiles:.c=.o}
SLAVEOBJS+=${ptrFiles:.c=.o}

libinterface='../lib/libinterface.a'

default:${libinterface}

${libinterface}:${OBJS} ${HOSTOBJS} ${SLAVEOBJS}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@
${OBJS}: %.o:%.cpp
	${CXX} ${CXXFALGS} ${INCLUDE} -c $< -o $@
${HOSTOBJS}: %.o:%.c 
	${CC} ${CHOSTFLAGS} ${INCLUDE} -c $< -o $@
${SLAVEOBJS}: %.o:%.c
	${CC} ${CSLAVEFLAGS} ${INCLUDE} -c $< -o $@

.PHONY:clean
clean:
	rm -f ${OBJS} ${HOSTOBJS} ${SLAVEOBJS} ${libinterface}

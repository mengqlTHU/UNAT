#############################################################
# path setting
#############################################################
PROJECT=/home/export/online1/systest/swrh/lhb/UNAT_Renhu
LIBPATH=${PROJECT}/lib
INCLUDE=${PROJECT}/include
OBJPATH=${PROJECT}/objects

#############################################################
# build tool set setting
#############################################################
CC= 
CXX=mpiCC 
SLAVECC= 
AR= 
RANLIB= 
LD= 


#############################################################
# source directory
#############################################################
src= test.cpp

##################################################################
# Targets
##################################################################
test: test.o ${LIBPATH}/libUNAT.a
	${CXX} -L${LIBPATH} -lUNAT  $< -o $@

test.o: test.cpp 
	${CXX} -I${INCLUDE} -c $< 

.PHONY:clean
clean:
	@rm -f test *.o


#############################################################
# path setting
#############################################################
PROJECT= 
LIBPATH= 
INCLUDE= 
OBJPATH= 

#############################################################
# build tool set setting
#############################################################
CC= 
CXX= 
SLAVECC= 
AR= 
RANLIB= 
LD= 

#############################################################
# source directory
#############################################################
EXEDIR:=./directSegment

#############################################################
# make
#############################################################
all:
	(cd ${path}; make CC=${CC} CXX=${CXX} OBJPATH=${OBJPATH} \
	INLCUDE=${INLCUDE} SLAVECC=${SLAVECC} RANLIB=${RANLIB} LD=${LD})


.PHONY:clean
clean:
	@for path in $(EXEDIR); do \
		(cd $$path; \
		make -f ../Makefile.rule clean); \
	done

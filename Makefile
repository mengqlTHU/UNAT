#############################################################
# path setting
#############################################################
PROJECT=${WM_PROJECT_DIR}/src/OSspecific/performance/UNAT
LIBPATH=${PROJECT}/lib
INCLUDE=${PROJECT}/include
OBJPATH=${PROJECT}/objects
LIBROOT=${WM_PROJECT_DIR}/platforms/linux64swgccDPInt32Opt/lib

#############################################################
# build tool set setting
#############################################################
CC="sw5cc -host"
CXX="mpiCC"
SLAVECC="sw5cc -slave -msimd"
AR=swar cru
RANLIB=swranlib
LD=${CXX}

#############################################################
# source directory
#############################################################
SRCPATH:=${PROJECT}/tools \
		 ${PROJECT}/topology \
		 ${PROJECT}/iterator \
		 ${PROJECT}/wrappedInterface

TESTPATH=${PROJECT}/test

#############################################################
# make
#############################################################
.PHONY:all libobjs makepath

${LIBPATH}/libUNAT.a: libUNAT.a
	@cp $< $@
	@echo -e "\033[32mBuild completed! \033[0m"
	@cp $< ${LIBROOT}

libUNAT.a: libobjs
	${AR} $@ ${OBJPATH}/*.o
	${RANLIB} $@ 
	
test: libUNAT.a
	(cd ${path}; make CC=${CC} CXX=${CXX} OBJPATH=${OBJPATH} \
	INLCUDE=${INLCUDE} SLAVECC=${SLAVECC} RANLIB=${RANLIB} LD=${LD})

libobjs: makepath
	@for path in $(SRCPATH); do \
		echo -e "\033[32mEntering path $$path \033[0m"; \
		(cd $$path; \
		make -f ../Makefile.rule CC=${CC} CXX=${CXX} OBJPATH=${OBJPATH} \
		INLCUDE=${INLCUDE} SLAVECC=${SLAVECC} RANLIB=${RANLIB}); \
	done

makepath:
	@if [ ! -d ${OBJPATH} ]; then \
		echo -e "\033[32mMaking directory ${OBJPATH} \033[0m"; \
		mkdir ${OBJPATH}; \
	fi
	@if [ ! -d ${INCLUDE} ]; then \
		echo -e "\033[32mMaking directory ${INCLUDE} \033[0m"; \
		mkdir ${INCLUDE}; \
	fi
	@if [ ! -d ${LIBPATH} ]; then \
		echo -e "\033[32mMaking directory ${LIBPATH} \033[0m"; \
		mkdir ${LIBPATH}; \
	fi

.PHONY:clean
clean:
	rm -rf libUNAT.a ${OBJPATH} ${INCLUDE} ${LIBPATH}
	@for path in $(SRCPATH); do \
		(cd $$path; \
		make -f ../Makefile.rule clean); \
	done

CC=sw5cc
CXX=swg++453
LINKER=swld453
CFLAGS=-host -fPIC -O3
gFLAGS=-slave -fPIC -msimd -g
CXXFLAGS=-fPIC -O3
SFLAGS=-slave -fPIC -msimd -O3
LINKFLAGS=-hybrid
PROJDIR=/home/export/online1/systest/swrh/lhb/UNAT_new

files= \
${PROJDIR}/iterator/multiLevelBlockIterator/edge2VertexIter/edge2VertexIter_host.o \
${PROJDIR}/iterator/multiLevelBlockIterator/vertex2EdgeIter/vertex2EdgeIter_host.o \
${PROJDIR}/iterator/multiLevelBlockIterator/edge2VertexIter/funcPointer.o \
${PROJDIR}/iterator/multiLevelBlockIterator/edge2VertexIter/userFunc_host.o \
${PROJDIR}/iterator/multiLevelBlockIterator/vertex2EdgeIter/userFuncUnsymm_host.o \
${PROJDIR}/iterator/multiLevelBlockIterator/BlockOrdering/BlockOrdering.o \
${PROJDIR}/iterator/multiLevelBlockIterator/BlockOrdering/BOrderUtils.o \
${PROJDIR}/iterator/multiLevelBlockIterator/extensibleArray/extensibleLabelArray.o \
${PROJDIR}/iterator/multiLevelBlockIterator/extensibleArray/extensibleScalarArray.o \

slaveFiles= \
${PROJDIR}/iterator/multiLevelBlockIterator/edge2VertexIter/edge2VertexIter_slave.o \
${PROJDIR}/iterator/multiLevelBlockIterator/edge2VertexIter/userFunc_slave.o \
${PROJDIR}/iterator/multiLevelBlockIterator/vertex2EdgeIter/vertex2EdgeIter_slave.o \
${PROJDIR}/iterator/multiLevelBlockIterator/vertex2EdgeIter/userFuncUnsymm_slave.o \

cppFiles= \
${PROJDIR}/test/multiLevelBlock/test.o 

regFiles= \
${PROJDIR}/RL_MPI/register.o

tools= \
${PROJDIR}/tools/slaveUtils.o

INCLUDE=-I${PROJDIR}/topology/ -I${PROJDIR}/iterator/ -I${PROJDIR}/iterator/multiLevelBlockIterator/ -I${PROJDIR}/ -I${PROJDIR}/RL_MPI/ -I${PROJDIR}/iterator/multiLevelBlockIterator/edge2VertexIter/ -I${PROJDIR}/tools/ -I${PROJDIR}/iterator/multiLevelBlockIterator/vertex2EdgeIter/
LIBRARY=-L${PROJDIR}/lib/LIB_SW

OBJECT=a.out

default: ${cppFiles} ${files} ${slaveFiles} ${regFiles} ${tools}
	${LINKER} ${LINKFLAGS} $^ -o ${OBJECT} -lm -lstdc++ ${LIBRARY} -lmetis -lrlmpi
${cppFiles}:%.o:%.cpp
	${CXX} ${CXXFLAGS} -c $^ -o $@ ${INCLUDE}
${files}:%.o:%.c
	${CC} ${CFLAGS} -c $< -o $@ ${INCLUDE}
${slaveFiles}:%.o:%.c
	${CC} ${SFLAGS} -c $< -o $@ ${INCLUDE}
${regFiles}:%.o:%.C
	${CC} -host -O3 -c $< -o $@ ${INCLUDE}
${tools}:%.o:%.c
	${CC} ${gFLAGS} -c $< -o $@ 
clean:
	rm -rf ${OBJECT} 
	find ${PROJDIR} -name "*.o" | xargs rm -rf

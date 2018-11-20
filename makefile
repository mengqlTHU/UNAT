CC=sw5cc
CXX=swg++453
LINKER=swld453
CFLAGS=-host -fPIC -O3
gFLAGS=-slave -fPIC -g
CXXFLAGS=-fPIC -O3
SFLAGS=-slave -fPIC -msimd -O3
LINKFLAGS=-hybrid

files= \
./iterator/multiLevelBlockIterator/edge2VertexIter_host.o \
./iterator/multiLevelBlockIterator/BlockOrdering/BlockOrdering.o \
./iterator/multiLevelBlockIterator/BlockOrdering/BOrderUtils.o \
./iterator/multiLevelBlockIterator/extensibleArray/extensibleLabelArray.o \
./iterator/multiLevelBlockIterator/extensibleArray/extensibleScalarArray.o 

slaveFiles= \
./iterator/multiLevelBlockIterator/edge2VertexIter_slave.o

cppFiles= \
./test.o \
./iterator/multiLevelBlockIterator/edge2VertexIter_reg.o 

INCLUDE=-I./topology/ -I./iterator/ -I./iterator/multiLevelBlockIterator/ -I./ -I./RL_MPI/
LIBRARY=-L./lib/LIB_SW

OBJECT=a.out

default: ${cppFiles} ${files} ${slaveFiles}
	${LINKER} ${LINKFLAGS} $^ -o ${OBJECT} -lm -lstdc++ ${LIBRARY} -lmetis -lrlmpi
${cppFiles}:%.o:%.cpp
	${CXX} ${CXXFLAGS} -c $^ -o $@ ${INCLUDE}
${files}:%.o:%.c
	${CC} ${CFLAGS} -c $< -o $@ ${INCLUDE}
${slaveFiles}:%.o:%.c
	${CC} ${SFLAGS} -c $< -o $@ ${INCLUDE}
clean:
	rm -rf ${OBJECT} 
	find . -name "*.o" | xargs rm -rf

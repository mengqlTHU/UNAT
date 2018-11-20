CC=gcc
CXX=g++
LINKER=g++
CFLAGS=-fPIC -O3
gFLAGS=-fPIC -g
CXXFLAGS=-fPIC -O3
SFLAGS=-O3

files= \
./iterator/multiLevelBlockIterator/BlockOrdering/BlockOrdering.o \
./iterator/multiLevelBlockIterator/BlockOrdering/BOrderUtils.o \
./iterator/multiLevelBlockIterator/extensibleArray/extensibleLabelArray.o \
./iterator/multiLevelBlockIterator/extensibleArray/extensibleScalarArray.o \

INCLUDE=-I./topology/ -I./iterator/ -I./iterator/multiLevelBlockIterator/ -I./
LIBRARY=-L./lib/LIB_GCC

OBJECT=a.out

default: test.o ${files}
	${LINKER} $^ -o ${OBJECT} -lm ${LIBRARY} -lmetis
test.o:%.o:%.cpp
	${CXX} ${CXXFLAGS} -c $^ -o $@ ${INCLUDE}
${files}:%.o:%.c
	${CC} ${CFLAGS} -c $< -o $@
clean:
	rm -rf ${OBJECT} 
	find . -name "*.o" | xargs rm -rf

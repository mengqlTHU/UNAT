#ifndef VERTEX2EDGE_HOST_H
#define VERTEX2EDGE_HOST_H

#include <stdlib.h>
#include "swMacro.h"
#include "iterator_struct.h"
#include "edge2VertexIter_slave.h"

#ifdef __cplusplus
extern "C"{
#endif

void vertex2EdgeIteration_host(Arrays* neighbourData, Arrays* vertexData,
			void(*operatorFunPointer_host) (MLBFunParameters *MLBFunParas),
			void(*operatorFunPointer_slave) (MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas);

#ifdef __cplusplus
}
#endif

#endif

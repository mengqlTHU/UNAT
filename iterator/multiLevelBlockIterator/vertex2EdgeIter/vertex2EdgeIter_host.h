#ifndef VERTEX2EDGE_HOST_H
#define VERTEX2EDGE_HOST_H

#include <stdlib.h>
#include "swMacro.h"
#include "iterator_struct.h"
#include "edge2VertexIter_slave.h"

#ifdef __cplusplus
extern "C"{
#endif

typedef struct
{
	MLBParameters* MLBParas;
	swFloat*       data;
	swFloat*       x;
	swFloat*       b;
	swFloat*       diag;
	swInt*         firstEdgeVertices;
	swInt*         vertexNeighbor;
	swInt          isXExist;
	swInt          spIndex;  
}v2EParameters;

void vertex2EdgeIteration_host(Arrays* neighbourData, Arrays* vertexData,
//			void(*operatorFunPointer_host) (MLBFunParameters *MLBFunParas),
//			void(*operatorFunPointer_slave) (MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas);

#ifdef __cplusplus
}
#endif

#endif

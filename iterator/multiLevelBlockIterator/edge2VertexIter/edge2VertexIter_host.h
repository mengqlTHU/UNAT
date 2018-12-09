#ifndef EDGE2VERTEXITER_HOST_H
#define EDGE2VERTEXITER_HOST_H

#include <stdlib.h>
#include "swMacro.h"
#include "iterator_struct.h"
#include "edge2VertexIter_slave.h"

typedef struct
{
	MLBParameters* MLBParas;
	swFloat*       lower;
	swFloat*       upper;
	swFloat*       x;
	swFloat*       b;
	swFloat*       diag;
	swInt*         owner;
	swInt*         neighbor;
	swInt          spIndex;
	swInt          isXExist;
	swInt          dimension;
	swInt          vertexNumber;
	swInt          edgeNumber;
}e2VParameters;

void edge2VertexIteration_host(Arrays* edgeData, Arrays* vertexData,
//			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
//			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas);
void edge2VertexIteration_init(Arrays* edgeData, Arrays* vertexData,
//			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
//			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas);

#endif

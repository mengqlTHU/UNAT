#ifndef EDGE2VERTEXITER_HOST_H
#define EDGE2VERTEXITER_HOST_H

#include <stdlib.h>
#include "swMacro.h"
#include "iterator_struct.h"
#include "edge2VertexIter_slave.h"


swInt spIndex, maxXNum, maxCell, maxEdge, cpeBlockNum, mshBlockNum;
swInt *blockStarts, *cellStarts;
swInt *ownNeiSendList, *owner, *neighbor;
swFloat *upper, *lower, *diag, *x, *b;
bool isXExist;
void (*operatorFunPointer_h)(MLBFunParameters *MLBFunParas);
void (*operatorFunPointer_s)(MLBFunParameters *MLBFunParas);

void edge2VertexIteration_host(Arrays* edgeData, Arrays* vertexData,
			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas);
void initOwnNeiSendList();
void edge2VertexIteration_init(Arrays* edgeData, Arrays* vertexData,
			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas);

#endif

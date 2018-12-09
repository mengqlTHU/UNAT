#ifndef ITERATOR_STRUCT_H
#define ITERATOR_STRUCT_H

#include "swMacro.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C"{
#endif

typedef struct
{
	swFloat* A1Ptr;
	swFloat* A2Ptr;
	swFloat* A3Ptr;
	swFloat* A4Ptr;
	swInt    num;
	swInt    dimension;
}Arrays;

typedef struct
{
	swInt* sOwner;
	swInt* rOwner;
	swInt* sNeighbor;
	swInt* rNeighbor;
	swInt* diagOwner;
	swInt* diagNeighbor;
}topoArrays;

typedef struct
{
	Arrays*     edgeData;
	Arrays*     vertexData;
	topoArrays* tArrays;
	swInt       count;
	swInt       k1;
	swInt       k2;
	swInt       k3;
	swInt       flag;
}MLBFunParameters;

typedef struct
{
	swInt* blockStarts;
	swInt* blockStartsUnsymm;
	swInt* vertexStarts;
	swInt* owner;
	swInt* neighbor;
	swInt* firstEdgeVertices;
	swInt* vertexNeighbor;
	swInt  cpeBlockNum;
	swInt  mshBlockNum;
	swInt  mtxBlockNum;
	swInt  maxXNum;
	swInt  maxCells;
	swInt  maxEdges;
	void   (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas);
	void   (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas);
}MLBParameters;

//typedef struct
//{
//	Arrays* edgeData;
//	Arrays* vertexData;
//	swInt*  owner;
//	swInt*  neighbor;
//	swInt   count;
//	swInt   k1;
//	swInt   k2;
//	swInt   flag;
//}MLBFunParameters;

#ifdef __cplusplus
}
#endif

#endif

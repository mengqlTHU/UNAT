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
}Arrays;

typedef struct
{
	swInt* blockStarts;
	swInt* vertexStarts;
	swInt* owner;
	swInt* neighbor;
	swInt  cpeBlockNum;
	swInt  mshBlockNum;
	swInt  mtxBlockNum;
	swInt  maxXNum;
	swInt  maxCells;
	swInt  maxEdges;
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

#ifndef EDGE2VERTEXITER_SLAVE_H
#define EDGE2VERTEXITER_SLAVE_H

#include <stdlib.h>
#include "iterator_struct.h"
#ifdef __cplusplus
extern "C"{
#endif

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
	swInt       flag;
}MLBFunParameters;

#ifdef __cplusplus
}
#endif

#endif

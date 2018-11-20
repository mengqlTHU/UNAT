#ifndef ITERATOR_STRUCT_H
#define ITERATOR_STRUCT_H

#include "swMacro.h"
#include <stdlib.h>

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

#endif

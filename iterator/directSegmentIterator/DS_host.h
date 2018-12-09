#ifndef DS_HOST_H
#define DS_HOST_H

#include "swMacro.h"
#include "iterator_struct.h"

#ifdef __cplusplus
extern "C"{
#endif

typedef struct
{
	swInt* segStarts;
	swInt* edgeStarts;
	swInt* owner;
	swInt* neighbor;
	swInt  segNum;
	swInt  subSegNum;
	swInt  maxEdges;
	swInt  vertexNumber;
	swInt  edgeNumber;
	void   (*operatorFunPointer_host)();
	void   (*operatorFunPointer_slave)();
}DSParameters;

typedef struct
{
	DSParameters* DSParas;
	Arrays*       upper;
	Arrays*       lower;
	Arrays*       diag;
	Arrays*       vertex;
	swInt         spIndex;
	swInt         isXExist;
	swInt         dimension;
}DS_e2VParameters;


#ifdef __cplusplus
}
#endif

#endif

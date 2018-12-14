#ifndef DIRECTSEGMENTITERATOR_H
#define DIRECTSEGMENTITERATOR_H

#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct 
{
	// topology relevent
	swInt* startVertex;
	swInt* endVertex;
	swInt segNum;
	swInt subSegNum;
	swInt *segStarts;
	swInt *edgeStarts;
	swInt *startSeg;
	swInt *endSeg;
	swInt  connectNum;
	swInt* segEdgeNum;
	swInt* edgeNeiSeg;
	swInt  vertexNumber;

	// Run-time data
	swInt spIndex;

	// computing data
	Arrays* backEdgeData;
	Arrays* frontEdgeData;
	Arrays* selfConnData;
	Arrays* vertexData;

	// computing pointers
	e2v_hostFunPtr fun_host;
	e2v_slaveFunPtr fun_slave;

} DS_edge2VertexPara;

#ifdef __cplusplus
}
#endif



#endif // DIRECTSEGMENTITERATOR_H

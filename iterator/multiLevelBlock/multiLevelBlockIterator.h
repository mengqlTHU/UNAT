#ifndef MULTILEVELBLOCKITERATOR_H
#define MULTILEVELBLOCKITERATOR_H

#include "iterator.h"

#ifdef __cplusplus
extern "C"
{
#include "RlmpiInitializer.hxx"
#endif

typedef struct 
{
	// topology relevent
	swInt* owner;
	swInt* neighbor;
	swInt* cellStarts;
	swInt* blockStarts;
	swInt* recvStarts;
	swInt  vertexNumber;
	swInt  edgeNumber;
	swInt  mshBlockNum;
	swInt  maxEdges;
	swInt  maxCells;

	// Run-time data
	swInt spIndex;

	// Register communication topology
	RlmpiInfo *schedule_data;

	// computing data
	Arrays* backEdgeData;
	Arrays* frontEdgeData;
	Arrays* selfConnData;
	Arrays* vertexData;

	// computing pointers
	e2v_hostFunPtr fun_host;
	e2v_slaveFunPtr fun_slave;

} MLB_edge2VertexPara;

#ifdef __cplusplus
}
#endif



#endif // MULTILEVELBLOCKITERATOR_H

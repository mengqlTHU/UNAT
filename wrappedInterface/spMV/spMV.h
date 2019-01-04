#ifndef SPMV_H
#define SPMV_H
#include <stdlib.h>
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
void slave_spMV_slave(Arrays* backEdgeData, Arrays* frontEdgeData,
			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
			swInt* endVertices);
void spMV_host(Arrays* backEdgeData, Arrays* frontEdgeData,
			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
			swInt* endVertices);
#ifdef __cplusplus
}
#endif

#endif

#ifndef INTEGRATE_H
#define INTEGRATE_H
#include <stdlib.h>
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
define_e2v_FunPtr(integrate);
define_e2v_slaveFunPtr(integrate);
//void slave_integrate(Arrays* backEdgeData, Arrays* frontEdgeData,
//			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
//			swInt* endVertices);
//void integrate(Arrays* backEdgeData, Arrays* frontEdgeData,
//			Arrays* selfConnData, Arrays* vertexData, swInt* startVertices,
//			swInt* endVertices);
#ifdef __cplusplus
}
#endif

#endif

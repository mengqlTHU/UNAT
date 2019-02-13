#ifndef INTERPOLATION_H
#define INTERPOLATION_H
#include <stdlib.h>
#include "swMacro.h"
#include "iterator.h"

//using namespace UNAT;
#ifdef __cplusplus
extern "C"
{
#endif
	//define 2 function pointers
define_e2v_FunPtr(interpolation);
define_e2v_slaveFunPtr(interpolation);
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

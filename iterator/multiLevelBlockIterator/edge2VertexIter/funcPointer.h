#ifndef FUNCPOINTER_H
#define FUNCPOINTER_H

#ifdef __cplusplus
extern "C"{
#endif

#include "edge2VertexIter_slave.h"

typedef void (*PF)(MLBFunParameters *MLBFunParas);
PF funcPointer_host();
PF funcPointer_slave();

#ifdef __cplusplus
}
#endif

#endif

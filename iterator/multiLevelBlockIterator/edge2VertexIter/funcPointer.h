#ifndef FUNCPOINTER_H
#define FUNCPOINTER_H

#ifdef __cplusplus
extern "C"{
#endif

#include "edge2VertexIter_slave.h"

typedef void (*PF)(MLBFunParameters *MLBFunParas);
PF funcPointer_host(swInt flag);
PF funcPointer_slave(swInt flag);

#ifdef __cplusplus
}
#endif

#endif

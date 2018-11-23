#ifndef USERFUNC_HOST_H
#define USERFUNC_HOST_H

#include <stdlib.h>
#include "swMacro.h"
#include "edge2VertexIter_slave.h"

#ifdef __cplusplus
extern "C"{
#endif

void spMV(MLBFunParameters *MLBFunParas);
void integrate(MLBFunParameters *MLBFunParas);

#ifdef __cplusplus
}
#endif

#endif

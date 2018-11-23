#ifndef USERFUNC_SLAVE_H
#define USERFUNC_SLAVE_H

#include <stdlib.h>
#include "swMacro.h"
#include "edge2VertexIter_slave.h"

#ifdef __cplusplus
extern "C"{
#endif

#define SLAVE_FUNC(funcname) slave_##funcname

void spMV(MLBFunParameters *MLBFunParas);
void integrate(MLBFunParameters *MLBFunParas);

#ifdef __cplusplus
}
#endif

#endif

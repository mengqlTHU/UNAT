#ifndef USERFUNC_SLAVE_H
#define USERFUNC_SLAVE_H

#include <stdlib.h>
#include "swMacro.h"
#include "iterator_struct.h"

#ifdef __cplusplus
extern "C"{
#endif

#define SLAVE_FUNC(funcname) slave_##funcname

void spMVUnsymm(MLBFunParameters *MLBFunParas);
void integrateUnsymm(MLBFunParameters *MLBFunParas);

#ifdef __cplusplus
}
#endif

#endif

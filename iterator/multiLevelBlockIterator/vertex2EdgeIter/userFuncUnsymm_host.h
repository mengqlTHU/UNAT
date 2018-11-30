#ifndef USERFUNCUNSYMM_HOST_H
#define USERFUNCUNSYMM_HOST_H

#include <stdlib.h>
#include "swMacro.h"
#include "iterator_struct.h"

#ifdef __cplusplus
extern "C"{
#endif

void spMVUnsymm(MLBFunParameters *MLBFunParas);
void integrateUnsymm(MLBFunParameters *MLBFunParas);

#ifdef __cplusplus
}
#endif

#endif

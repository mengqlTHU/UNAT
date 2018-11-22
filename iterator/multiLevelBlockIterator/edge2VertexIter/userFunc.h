#ifndef USERFUNC_H
#define USERFUNC_H

#include <stdlib.h>
#include "swMacro.h"
#include "iterator_struct.h"

#ifdef __cplusplus
extern "C"{
#endif

#define SLAVE_FUNC(funcname) slave_##funcname

void spMV(MLBFunParameters *MLBFunParas);

#ifdef __cplusplus
}
#endif

#endif

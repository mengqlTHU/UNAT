#ifndef FUNCPOINTER_H
#define FUNCPOINTER_H

#ifdef __cplusplus
extern "C"{
#endif

#include "iterator_struct.h"

typedef void (*PF)(MLBFunParameters *MLBFunParas);
PF funcPointer();

#ifdef __cplusplus
}
#endif

#endif

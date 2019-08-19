#ifndef SW_MLBINTERFACE_H
#define SW_MLBINTERFACE_H
#include "swMacro.h"

extern "C" void constructMLBIterator(swInt* rowAddr, swInt* colAddr, 
			swInt* cellWeights, swInt* faceWeights, swInt faceNum);
extern "C" swInt* getBlockStarts();
extern "C" swInt* getPostEdgeOrder();
extern "C" swInt* getPostVertexOrder();

#endif

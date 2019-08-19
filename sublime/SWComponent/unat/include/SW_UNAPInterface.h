#ifndef SW_UNAPINTERFACE_H
#define SW_UNAPINTERFACE_H

#include "multiLevelBlockIterator.hpp"

void constructUNAPIterator(MultilevelBlockIterator *mlbIter,
			swInt* rowAddr, swInt* colAddr, 
			swInt* cellWeights, swInt* faceWeights, swInt faceNum);
swInt* getPostEdgeOrder();
swInt* getPostVertexOrder();

#endif

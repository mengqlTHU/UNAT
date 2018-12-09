#ifndef DS_SLAVE_H
#define DS_SLAVE_H

#include "swMacro.h"

// Topology
__thread_local_fix swInt   *segStarts,*edgeStarts,*owner,*neighbor;
// Field data
__thread_local_fix Arrays *upper,*lower,*diag,*vertex;
// spIndex: index of segment
// maxEdges: maximum count of edges among the subsegments in same segment
// isXExist: whether x/b is contained in vertex data.
// dimension: the dimension of data; 1--x; 2--x,y; 3--x,y,z;
__thread_local_fix swInt   spIndex,maxEdges,cpeBlockNum,isXExist,dimension;
//
__thread_local_fix swInt   vertexNumber,edgeNumber;

// variables(***_slave) are allocated in LDM
// sVertex_slave: local vertex data;
// rVertex_slave: global vertex data including received data from other
// blocks for the arrays which fArratInOut[i]=COPYIN, while others keep
// local.
__thread_local_fix Arrays  *sVertex_slave,*rVertex_slave,*diag_slave;
__thread_local_fix swInt   *owner_slave,*neighbor_slave;
__thread_local_fix swInt   *segStarts_slave,*edgeStarts_slave;
// colIdx: the column block index of each edge
// colNum: the counts of each column block
// accuColNum: another representation of colNum
__thread_local_fix swInt   *colIdx,*colNum,*accuColNum,*colGloablIdx;
// sLower_slave: lower data stored upper triangle.
// rLower_slave: lower data received from other blocks
__thread_local_fix Arrays  *sLower_slave,*upper_slave,*rLower_slave;
// record the base in quicksort
__thread_local_fix swInt   pivOwner,pivNeighbor;
__thread_local_fix swFloat **pivLower,**pivUpper;

#endif

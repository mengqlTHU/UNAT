#include <stdlib.h>
#include <stdio.h>
#include "slave.h"
#include "dma.h"
#include "rlmpi.h"
#include "swMacro.h"
#include "iterator_struct.h"
#include "userFunc_slave.h"

__thread_local volatile int myId;

void vertex2EdgeIter_slave()
{
	myId = athread_get_id(-1);
	if(myId==0) printf("hello world!\n");

	swFloat sendX_slave[maxXNum],recvXNum[maxXNum],b_slave[maxXNum];
	swFloat diag[maxCells*2],nonDiag[maxEdges*2];
	swInt cell_slave[BLOCKNUM64K+1];
	swInt blockStarts_slave[BLOCKNUM64K*4];
}

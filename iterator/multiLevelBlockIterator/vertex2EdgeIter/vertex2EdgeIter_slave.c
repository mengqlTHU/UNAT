#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "slave.h"
#include "dma.h"
#include "rlmpi.h"
#include "swMacro.h"
#include "iterator_struct.h"
#include "userFunc_slave.h"
#include "vertex2EdgeIter_host.h"

__thread_local volatile int myId;

//extern void (*operatorFunPointer_s)(MLBFunParameters *MLBFunParas);
//extern swInt spIndex,maxXNum,maxCells,maxEdges,cpeBlockNum;
//extern swInt *blockStarts,*cellStarts,*ownNeiSendList;
//extern swFloat *data,*x,*b,*diag;

void vertex2EdgeIter_slave(v2EParameters *v2EParas)
{
	myId = athread_get_id(-1);

	swInt maxXNum,maxCells,maxEdges,cpeBlockNum,spIndex;
	spIndex     = v2EParas->spIndex;
	maxXNum     = v2EParas->MLBParas->maxXNum;
	maxCells    = v2EParas->MLBParas->maxCells;
	maxEdges    = v2EParas->MLBParas->maxEdges;
	cpeBlockNum = v2EParas->MLBParas->cpeBlockNum;

	swInt *blockStarts,*cellStarts,*firstEdgeVertices,*vertexNeighbor;
	blockStarts       = v2EParas->MLBParas->blockStartsUnsymm;
	cellStarts        = v2EParas->MLBParas->vertexStarts;
	firstEdgeVertices = v2EParas->MLBParas->firstEdgeVertices;
	vertexNeighbor    = v2EParas->MLBParas->vertexNeighbor;

	swFloat *data,*x,*b,*diag;
	data = v2EParas->data;
	x    = v2EParas->x;
	b    = v2EParas->b;
	diag = v2EParas->diag;

//	swInt totalLength
//		= maxXNum*2*sizeof(swFloat)
//		+ maxEdges*2*sizeof(swFloat)
//		+ maxEdges*2*sizeof(swInt)
//		+ BLOCKNUM64K*5*sizeof(swInt);
//	if(totalLength>64*1024*0.94 && myId==0)
//	{
//		printf("The allocated LDM exceeds 64KB, the memory size is %d\n",
//					totalLength);
//	}
	swFloat sendX_slave[maxXNum],b_slave[maxXNum];
	swFloat data_slave[maxEdges],recvX_slave[maxEdges];
	swInt firstEdgeVertices_slave[maxEdges];
	swInt vertexNeighbor_slave[maxEdges];
	swInt cellStart,cellEnd;
	swInt blockStarts_slave[BLOCKNUM64K*4],ownNeiSendIdx[BLOCKNUM64K];

	swInt startIdx,startBlockIdx,blockLen,endIdx;
	swInt i,j,recvLower,recvUpper;
	swInt row = myId+spIndex*BLOCKNUM64K;
	startBlockIdx = row*cpeBlockNum+spIndex*BLOCKNUM64K;
	cellStart = cellStarts[row];
	cellEnd   = cellStarts[row+1];
	DMA_Get(&blockStarts_slave[0],&blockStarts[4*startBlockIdx],
				BLOCKNUM64K*4*sizeof(swInt));
	DMA_Get(&b_slave[0],&b[cellStart],(cellEnd-cellStart)*sizeof(swFloat));

	blockLen = blockStarts_slave[4*(BLOCKNUM64K-1)+3]-blockStarts_slave[2];
	DMA_Get(&data_slave[0],&data[blockStarts_slave[2]],
				blockLen*sizeof(swFloat));
	DMA_Get(&firstEdgeVertices_slave[0],
				&firstEdgeVertices[blockStarts_slave[2]],
				blockLen*sizeof(swInt));

	Arrays neighborData_slave = {&data_slave[0],NULL,NULL,NULL,0};
	Arrays vertexData_slave = {&b_slave[0],&sendX_slave[0],
		&recvX_slave[0],NULL,0};
	topoArrays tArrays = {&firstEdgeVertices_slave[0],NULL,
		&vertexNeighbor_slave[0],NULL,NULL,NULL};
	MLBFunParameters MLBFunParas = {&neighborData_slave, &vertexData_slave,
		&tArrays, 0, 0, 0, 0};

	if(x!=NULL)
	{
	DMA_Get(&sendX_slave[0],&x[cellStart],
				(cellEnd-cellStart)*sizeof(swFloat));
	DMA_Get(&vertexNeighbor_slave[0], &vertexNeighbor[blockStarts_slave[2]],
				blockLen*sizeof(swInt));
	DMA_Get(&recvX_slave[0],&diag[cellStart],
				(cellEnd-cellStart)*sizeof(swFloat));
	MLBFunParas.count = cellEnd-cellStart;
	MLBFunParas.flag  = 1;
	v2EParas->MLBParas->operatorFunPointer_slave(&MLBFunParas);

	for(i=0;i<cellEnd-cellStart;i++)
	{
//		if(i+cellStart==10723) printf("slave diag:%d,%f,%f\n",i,recvX_slave[i],sendX_slave[i]);
		b_slave[i] += recvX_slave[i]*sendX_slave[i];
	}


	for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
	for(i=0;i<_total_send_pcg;i++)
	{
		startIdx = blockStarts_slave[4*_sPacks[i].dst_id+2]
			+ ownNeiSendIdx[_sPacks[i].dst_id]-blockStarts_slave[2];
		endIdx   = blockStarts_slave[4*_sPacks[i].dst_id+3]
			- blockStarts_slave[2];
		for(j=0;j<6;j++)
		{
			if(startIdx+j>=endIdx) break;
			_sPacks[i].data[j] =
				sendX_slave[firstEdgeVertices_slave[startIdx+j]-cellStart];
			_sPacks[i].indM = j+1;
		}
		ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].indM;
	}
	transform_data();
	recvLower = 0;
	recvUpper = 0;
	for(i=0;i<_total_recv_pcg;i++)
	{
		if(_rPacks[i].src_id>myId)
		{
			for(j=0;j<_rPacks[i].indM;j++)
			{
				recvX_slave[recvUpper+blockStarts_slave[4*myId+3]
					- blockStarts_slave[2]] = _rPacks[i].data[j];
				recvUpper++;
			}
		}else if(_rPacks[i].src_id<myId)
		{
			for(j=0;j<_rPacks[i].indM;j++)
			{
				recvX_slave[recvLower] = _rPacks[i].data[j];
				recvLower++;
			}
		}
	}
	for(i=blockStarts_slave[4*myId+2];i<blockStarts_slave[4*myId+3];i++)
	{
		recvX_slave[i-blockStarts_slave[2]]
			= sendX_slave[vertexNeighbor_slave[i-blockStarts_slave[2]]
			- cellStart];
	}
	}

	MLBFunParas.count = blockLen;
	MLBFunParas.k1    = cellStart;
	v2EParas->MLBParas->operatorFunPointer_slave(&MLBFunParas);

//	for(i=0;i<blockLen;i++)
//	{
////		if(firstEdgeVertices_slave[i]==10723) printf("slave nondiag:%d,%f,%f\n",i,data_slave[i],recvX_slave[i]);
////		if(firstEdgeVertices_slave[i]>=cellEnd) printf("%d,%d,%d\n",i,firstEdgeVertices_slave[i],cellEnd);
////		assert(firstEdgeVertices_slave[i]<cellEnd);
//		b_slave[firstEdgeVertices_slave[i]-cellStart]
//			+= data_slave[i]*recvX_slave[i];
//	}


	DMA_Put(&b[cellStart],&b_slave[0],(cellEnd-cellStart)*sizeof(swFloat));
}

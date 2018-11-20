#include "edge2VertexIter_host.h"
#include <athread.h>
#include <stdlib.h>
#include <stdio.h>

extern SLAVE_FUN(func)();

void edge2VertexIteration_init(Arrays* edgeData, Arrays* vertexData,
			void (*operatorFunPointer)(Arrays*, Arrays*, swInt*, swInt*),
			MLBParameters* MLBParas)
{
	maxXNum     = MLBParas->maxXNum;
	maxCell     = MLBParas->maxCells;
	maxEdge     = MLBParas->maxEdges;
	cpeBlockNum = MLBParas->cpeBlockNum;
	mshBlockNum = MLBParas->mshBlockNum;
	blockStarts = MLBParas->blockStarts;
	cellStarts  = MLBParas->vertexStarts;
	owner       = MLBParas->owner;
	neighbor    = MLBParas->neighbor;
	lower       = edgeData->A1Ptr;
	upper       = edgeData->A2Ptr;
	b           = vertexData->A1Ptr;
	x           = vertexData->A2Ptr;
	diag        = vertexData->A3Ptr;
	operatorFunPointer_host = operatorFunPointer;

	initOwnNeiSendList();
	athread_init();
}

void edge2VertexIteration_host(Arrays* edgeData, Arrays* vertexData,
			void (*operatorFunPointer)(Arrays*, Arrays*, swInt*, swInt*),
			MLBParameters* MLBParas)
{
	edge2VertexIteration_init(edgeData, vertexData,
				operatorFunPointer, MLBParas);
	athread_spawn(func, 0);
	athread_join();
	athread_halt();
	printf("edge2VertexIteration_host\n");
}

void initOwnNeiSendList(){
	int num = mshBlockNum*BLOCKNUM64K*BLOCKNUM64K;
	int startIdx,endIdx,startListIdx,i,j,k,row,col,blockIdx,upperIdx,lowerIdx;
	ownNeiSendList = (int*)malloc(sizeof(int)*num);
	for(k=0;k<num;k++){ownNeiSendList[k]=0;}
	for(i=0;i<mshBlockNum;i++){
		for(j=0;j<BLOCKNUM64K;j++){
			row = i*BLOCKNUM64K+j;
			col = row+1;
			while(col < (i+1)*BLOCKNUM64K && col > i*BLOCKNUM64K){
				blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
				startIdx = blockStarts[4*blockIdx+2];
				endIdx   = blockStarts[4*blockIdx+3];
				upperIdx = row * BLOCKNUM64K + col - i * BLOCKNUM64K;
				lowerIdx = col * BLOCKNUM64K + row - i * BLOCKNUM64K;
				col++;
				if(endIdx-startIdx==0) continue;
				ownNeiSendList[upperIdx] = (endIdx - startIdx-1)/6+1;
				ownNeiSendList[lowerIdx] = (endIdx - startIdx-1)/6+1;
			}
		}
	}
	for(i=0;i<mshBlockNum;i++){
		for(j=0;j<BLOCKNUM64K;j++){
			row = i*BLOCKNUM64K+j;
			for(col=0;col<BLOCKNUM64K;col++){
				upperIdx = row * BLOCKNUM64K + col;
				lowerIdx = (col + (i * BLOCKNUM64K)) * BLOCKNUM64K + row - (i*BLOCKNUM64K);
				if(ownNeiSendList[upperIdx] != ownNeiSendList[lowerIdx]){
					printf("%d,%d,%d,%d\n",row,col,ownNeiSendList[upperIdx],ownNeiSendList[lowerIdx]);
				}
			}
		}
	}


}

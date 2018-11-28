#include "vertex2EdgeIter_host.h"
#include <sys/time.h>
#include <athread.h>
#include <stdio.h>

extern SLAVE_FUN(vertex2EdgeIter_slave)();
swInt spIndex, cpeBlockNum, mshBlockNum, mtxBlockNum,
	  maxXNum, maxCells, maxEdges;
swInt *blockStarts, *cellStarts, *ownNeiSendList;
swInt *firstEdgeVertices, *vertexNeighbor;
swFloat *data, *x, *b, *diag;
void initOwnNeiSendListV2E();
void (*operatorFunPointer_h)(MLBFunParameters *MLBFunParas);
void (*operatorFunPointer_s)(MLBFunParameters *MLBFunParas);

void vertex2EdgeIteration_host(Arrays* neighbourData, Arrays* vertexData,
			void(*operatorFunPointer_host) (MLBFunParameters *MLBFunParas),
			void(*operatorFunPointer_slave) (MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas)
{
	maxXNum              = MLBParas->maxXNum;
	maxCells             = MLBParas->maxCells;
	maxEdges             = MLBParas->maxEdges;
	cpeBlockNum          = MLBParas->cpeBlockNum;
	mtxBlockNum          = MLBParas->mtxBlockNum;
	mshBlockNum          = MLBParas->mshBlockNum;
	blockStarts          = MLBParas->blockStarts;
	cellStarts           = MLBParas->vertexStarts;
	data                 = neighbourData->A1Ptr;
	b                    = vertexData->A1Ptr;
	x                    = vertexData->A2Ptr;
	diag                 = vertexData->A3Ptr;
	firstEdgeVertices    = MLBParas->firstEdgeVertices;
	vertexNeighbor       = MLBParas->vertexNeighbor;
	operatorFunPointer_h = operatorFunPointer_host;
	operatorFunPointer_s = operatorFunPointer_slave;

	initOwnNeiSendListV2E();

	athread_init();
	athread_spawn(vertex2EdgeIter_slave,0);
	athread_join();
	athread_halt();


	int i,j,k,row,col,startIdx,endIdx,blockIdx;
	struct timeval start,end;
	gettimeofday(&start,NULL);

	for(i=0;i<vertexData->num;i++)
	{
		b[i] = diag[i]*x[i];
	}

	for(i=0;i<cpeBlockNum;i++)
	{
		for(j=0;j<cpeBlockNum;j++)
		{
			blockIdx = i*cpeBlockNum+j;
			startIdx = blockStarts[4*blockIdx+2];
			endIdx   = blockStarts[4*blockIdx+3];
			if(i==j)
			{
				for(k=startIdx;k<endIdx;k++)
				{
					b[firstEdgeVertices[k]]+=data[k]*x[vertexNeighbor[k]];
				}
			} else if(i<j)
			{
				for(k=startIdx;k<endIdx;k++)
				{
					b[firstEdgeVertices[k]]+=data[k]*x[vertexNeighbor[k]];
				}
			} else if(i>j)
			{
				for(k=startIdx;k<endIdx;k++)
				{
					b[firstEdgeVertices[k]]+=data[k]*x[vertexNeighbor[k]];
				}
			}
		}
	}

	gettimeofday(&end,NULL);
	int timeuse
		= 1000000*(end.tv_sec-start.tv_sec)+end.tv_usec-start.tv_usec;
	printf("reordered CPU Processor time: %f us\n",(double)timeuse);
	LOG("vertex2EdgeIteration_host");

}


void initOwnNeiSendListV2E()
{
	int num = mshBlockNum*BLOCKNUM64K*BLOCKNUM64K;
	int startIdx,endIdx,i,j,k,row,col,blockIdx,index;
	ownNeiSendList = (swInt*)malloc(sizeof(swInt)*num);
	for(i=0;i<num;i++) {ownNeiSendList[i]=0;}
	for(i=0;i<mshBlockNum;i++)
	{
		for(j=0;j<BLOCKNUM64K;j++)
		{
			row = i*BLOCKNUM64K+j;
			for(k=i*BLOCKNUM64K;k<(i+1)*BLOCKNUM64K;k++)
			{
				if(row==k) 	continue;
				blockIdx = row*cpeBlockNum+k;
				startIdx = blockStarts[4*blockIdx+2];
				endIdx   = blockStarts[4*blockIdx+3];
				index = row*BLOCKNUM64K+k-i*BLOCKNUM64K;
				if(endIdx-startIdx==0) continue;
				ownNeiSendList[index] = (endIdx-startIdx-1)/6+1;
			}
		}
	}
}

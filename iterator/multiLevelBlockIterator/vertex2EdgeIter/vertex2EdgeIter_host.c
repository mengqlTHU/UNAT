#include "vertex2EdgeIter_host.h"
#include <sys/time.h>
#include <athread.h>
#include <stdio.h>

extern SLAVE_FUN(vertex2EdgeIter_slave)(v2EParameters *v2EParas);
//swInt spIndex, cpeBlockNum, mshBlockNum, mtxBlockNum,
//	  maxXNum, maxCells, maxEdges;
//swInt *blockStarts, *cellStarts, *ownNeiSendList;
//swInt *firstEdgeVertices, *vertexNeighbor;
//swFloat *data, *x, *b, *diag;
void initOwnNeiSendListV2E(v2EParameters *v2EParas, swInt* ownNeiSendList);
//void (*operatorFunPointer_h)(MLBFunParameters *MLBFunParas);
//void (*operatorFunPointer_s)(MLBFunParameters *MLBFunParas);

void vertex2EdgeIteration_host(Arrays* neighbourData, Arrays* vertexData,
//			void(*operatorFunPointer_host) (MLBFunParameters *MLBFunParas),
//			void(*operatorFunPointer_slave) (MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas)
{
	swInt   cpeBlockNum,i,j,blockIdx,startIdx,endIdx,idxMP,spIndex,row,col;
	swInt   *firstEdgeVertices,*vertexNeighbor;
	swInt   *blockStarts,*blockStartsUnsymm,*owner,*neighbor,*cellStarts;
	swFloat *data,*b,*x,*diag;
//	maxXNum              = MLBParas->maxXNum;
//	maxCells             = MLBParas->maxCells;
//	maxEdges             = MLBParas->maxEdges;
	cpeBlockNum          = MLBParas->cpeBlockNum;
//	mtxBlockNum          = MLBParas->mtxBlockNum;
//	mshBlockNum          = MLBParas->mshBlockNum;
	blockStarts          = MLBParas->blockStarts;
	blockStartsUnsymm    = MLBParas->blockStartsUnsymm;
	cellStarts           = MLBParas->vertexStarts;
	data                 = neighbourData->A1Ptr;
	b                    = vertexData->A1Ptr;
	x                    = vertexData->A2Ptr;
	diag                 = vertexData->A3Ptr;
	firstEdgeVertices    = MLBParas->firstEdgeVertices;
	vertexNeighbor       = MLBParas->vertexNeighbor;
	owner                = MLBParas->owner;
	neighbor             = MLBParas->neighbor;
//	operatorFunPointer_h = operatorFunPointer_host;
//	operatorFunPointer_s = operatorFunPointer_slave;
	v2EParameters v2EParas;
	v2EParas.MLBParas = MLBParas;
	v2EParas.data     = neighbourData->A1Ptr;
	v2EParas.b        = vertexData->A1Ptr;
	v2EParas.x        = vertexData->A2Ptr;
	v2EParas.diag     = vertexData->A3Ptr;

	athread_init();

	swInt mshBlockNum = v2EParas.MLBParas->mshBlockNum;
	swInt num = mshBlockNum*BLOCKNUM64K*BLOCKNUM64K;
	swInt* ownNeiSendList = (swInt*)malloc(sizeof(swInt)*num);
	initOwnNeiSendListV2E(&v2EParas, ownNeiSendList);
	initSendList(ownNeiSendList,mshBlockNum);

	struct timeval start,end;
	gettimeofday(&start,NULL);
	for(i=0;i<BLOCKNUM64K;i++)
	{
		blockIdx = i*cpeBlockNum+BLOCKNUM64K;
		startIdx = blockStartsUnsymm[4*blockIdx+2];
		blockIdx = (i+1)*cpeBlockNum-1;
		endIdx   = blockStartsUnsymm[4*blockIdx+3];
		for(j=startIdx;j<endIdx;j++)
		{
			if(firstEdgeVertices[j]==10723) printf("1:%d,%f,%f\n",j,data[j],x[vertexNeighbor[j]]);
			b[firstEdgeVertices[j]] += data[j]*x[vertexNeighbor[j]];
		}
	}

	for(spIndex=0;spIndex<mshBlockNum;spIndex++)
	{
		initTable(spIndex);
		v2EParas.spIndex  = spIndex;
		athread_spawn(vertex2EdgeIter_slave,&v2EParas);
		idxMP = spIndex+1;
		if(idxMP<mshBlockNum)
		{
    		for(i=0;i<BLOCKNUM64K;i++)
    		{
    			row = idxMP*BLOCKNUM64K+i;
    			col = 0;
    			blockIdx = row*cpeBlockNum+col;
    			startIdx = blockStartsUnsymm[4*blockIdx+2];
    			col = idxMP*BLOCKNUM64K-1;
    			blockIdx = row*cpeBlockNum+col;
    			endIdx   = blockStartsUnsymm[4*blockIdx+3];
    			for(j=startIdx;j<endIdx;j++)
    			{
			if(firstEdgeVertices[j]==10723) printf("lower:%d,%f,%f\n",j,data[j],x[vertexNeighbor[j]]);
    				b[firstEdgeVertices[j]] += data[j]*x[vertexNeighbor[j]];
    			}
    
    			col = (idxMP+1)*BLOCKNUM64K;
    			if(col>=cpeBlockNum) continue;
    			blockIdx = row*cpeBlockNum+col;
    			startIdx = blockStartsUnsymm[4*blockIdx+2];
    			col = cpeBlockNum-1;
    			blockIdx = row*cpeBlockNum+col;
    			endIdx   = blockStartsUnsymm[4*blockIdx+2];
    			for(j=startIdx;j<endIdx;j++)
    			{
			if(firstEdgeVertices[j]==10723) printf("upper:%d,%f,%f\n",j,data[j],x[vertexNeighbor[j]]);
    				b[firstEdgeVertices[j]] += data[j]*x[vertexNeighbor[j]];
    			}
    		}
		}
		athread_join();
		destroyTable(spIndex);
	}
	gettimeofday(&end,NULL);
	int timeuse
		= 1000000*(end.tv_sec-start.tv_sec)+end.tv_usec-start.tv_usec;
	printf("v2e: Slave Processor time: %f us\n",(double)timeuse);

	athread_halt();


//	int i,j,k,row,col,startIdx,endIdx,blockIdx;
//	struct timeval start,end;
//	gettimeofday(&start,NULL);
//
//	for(i=0;i<vertexData->num;i++)
//	{
//		b[i] = diag[i]*x[i];
//	}
//
//	for(i=0;i<cpeBlockNum;i++)
//	{
//		for(j=0;j<cpeBlockNum;j++)
//		{
//			blockIdx = i*cpeBlockNum+j;
//			startIdx = blockStarts[4*blockIdx+2];
//			endIdx   = blockStarts[4*blockIdx+3];
//			if(i==j)
//			{
//				for(k=startIdx;k<endIdx;k++)
//				{
//					b[firstEdgeVertices[k]]+=data[k]*x[vertexNeighbor[k]];
//				}
//			} else if(i<j)
//			{
//				for(k=startIdx;k<endIdx;k++)
//				{
//					b[firstEdgeVertices[k]]+=data[k]*x[vertexNeighbor[k]];
//				}
//			} else if(i>j)
//			{
//				for(k=startIdx;k<endIdx;k++)
//				{
//					b[firstEdgeVertices[k]]+=data[k]*x[vertexNeighbor[k]];
//				}
//			}
//		}
//	}
//
//	gettimeofday(&end,NULL);
//	int timeuse
//		= 1000000*(end.tv_sec-start.tv_sec)+end.tv_usec-start.tv_usec;
//	printf("reordered CPU Processor time: %f us\n",(double)timeuse);
	LOG("vertex2EdgeIteration_host");

}


void initOwnNeiSendListV2E(v2EParameters *v2EParas, swInt* ownNeiSendList)
{
	swInt mshBlockNum = v2EParas->MLBParas->mshBlockNum;
	swInt cpeBlockNum = v2EParas->MLBParas->cpeBlockNum;
	swInt* blockStarts = v2EParas->MLBParas->blockStartsUnsymm;
	int num = mshBlockNum*BLOCKNUM64K*BLOCKNUM64K;
	int startIdx,endIdx,i,j,k,row,col,blockIdx,index;
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

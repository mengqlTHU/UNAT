#include "edge2VertexIter_host.h"
#include <athread.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

extern SLAVE_FUN(func)();

void edge2VertexIteration_init(Arrays* edgeData, Arrays* vertexData,
			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
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
	isXExist    = true;
	diag        = vertexData->A3Ptr;
	operatorFunPointer_h = operatorFunPointer_host;
	operatorFunPointer_s = operatorFunPointer_slave;
	if(x == NULL)
	{
		int i=0;
		isXExist=false;
		x = (swFloat*)malloc(vertexData->num*sizeof(swFloat));
		for(i=0;i<vertexData->num;i++){x[i]=1;}
	}

	initOwnNeiSendList();
	athread_init();
	initSendList(ownNeiSendList,mshBlockNum);
}

void edge2VertexIteration_host(Arrays* edgeData, Arrays* vertexData,
			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas)
{
	edge2VertexIteration_init(edgeData, vertexData, operatorFunPointer_host,
				operatorFunPointer_slave, MLBParas);

	swInt j,k,idxMP,blockIdx,startIdx,endIdx,row,col;
	topoArrays tArrays = {owner, NULL, neighbor, NULL, NULL, NULL};
	MLBFunParameters MLBFunParas = {edgeData, vertexData, &tArrays,
		0, 0, 0, 5};
	struct timeval start, end;
	gettimeofday(&start,NULL);
	for(j=0;j<BLOCKNUM64K;j++)
	{
		row = j;
		col = BLOCKNUM64K;
		blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
		startIdx = blockStarts[4*blockIdx+2];
		col = cpeBlockNum-1;
		blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
		endIdx = blockStarts[4*blockIdx+3];
		MLBFunParas.count = endIdx-startIdx;
		MLBFunParas.k1    = startIdx;
		operatorFunPointer_h(&MLBFunParas);
//		for(k=startIdx;k<endIdx;k++)
//		{
//			b[neighbor[k]] += lower[k] * x[owner[k]];
//			b[owner[k]]    += upper[k] * x[neighbor[k]];
//		}
	}

	for(spIndex=0;spIndex<mshBlockNum;spIndex++)
	{
		initTable(spIndex);
		athread_spawn(func, 0);
		idxMP = spIndex+1;
		for(j=0;j<BLOCKNUM64K;j++)
		{
			row      = idxMP*BLOCKNUM64K+j;
			col      = (idxMP+1)*BLOCKNUM64K;
//			col      = row;
			if(col>=cpeBlockNum) break;
			blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
			startIdx = blockStarts[4*blockIdx+2];
			col      = cpeBlockNum-1;
			blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
			endIdx   = blockStarts[4*blockIdx+3];
     		MLBFunParas.count = endIdx-startIdx;
	    	MLBFunParas.k1    = startIdx;
		    operatorFunPointer_h(&MLBFunParas);
//			for(k=startIdx;k<endIdx;k++)
//			{
//				b[neighbor[k]] +=lower[k] * x[owner[k]];
//				b[owner[k]]    +=upper[k] * x[neighbor[k]];
//			}
		}
		athread_join();
		destroyTable(spIndex);
	}
	gettimeofday(&end,NULL);
	int timeuse = 1000000*(end.tv_sec-start.tv_sec)
		+ end.tv_usec-start.tv_usec;
	printf("Slave Processor time: %f us\n",(double)timeuse);
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

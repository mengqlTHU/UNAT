#include "edge2VertexIter_host.h"
#include <athread.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

extern SLAVE_FUN(edge2VertexIter_slave)(e2VParameters *e2VParas);
//void (*operatorFunPointer_h)(MLBFunParameters *MLBFunParas);
//void (*operatorFunPointer_s)(MLBFunParameters *MLBFunParas);
void initOwnNeiSendList(MLBParameters *MLBParas, swInt* ownNeiSendList);

void edge2VertexIteration_init(Arrays* edgeData, Arrays* vertexData,
//			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
//			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas)
{
}

void edge2VertexIteration_host(Arrays* edgeData, Arrays* vertexData,
//			void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas),
//			void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas)
{
	swInt spIndex, maxXNum, maxCell, maxEdge, cpeBlockNum, mshBlockNum;
	swInt *blockStarts, *cellStarts, *owner, *neighbor, *ownNeiSendList;
	swFloat *upper, *lower, *diag, *x, *b;
	swInt isXExist, dimension, vertexNumber, edgeNumber;
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
	isXExist    = 1;
	diag        = vertexData->A3Ptr;
	dimension   = vertexData->dimension;
	vertexNumber= vertexData->num;
	edgeNumber  = edgeData->num;
//	operatorFunPointer_h = operatorFunPointer_host;
//	operatorFunPointer_s = operatorFunPointer_slave;
	if(x == NULL)
	{
		int i=0;
		isXExist=0;
		x = (swFloat*)malloc(vertexData->num*sizeof(swFloat)*dimension);
		for(i=0;i<vertexData->num*dimension;i++){x[i]=1;}
	}

	e2VParameters e2VParas = {MLBParas,lower,upper,x,b,diag,
		owner,neighbor,0,isXExist,dimension,vertexNumber,edgeNumber};

	int num = mshBlockNum*BLOCKNUM64K*BLOCKNUM64K;
	ownNeiSendList = (int*)malloc(sizeof(int)*num);
	initOwnNeiSendList(MLBParas, ownNeiSendList);
	athread_init();
	initSendList(ownNeiSendList,mshBlockNum);

	swInt j,k,idxMP,blockIdx,startIdx,endIdx,row,col;
	topoArrays tArrays = {owner, NULL, neighbor, NULL, NULL, NULL};
	MLBFunParameters MLBFunParas = {edgeData, vertexData, &tArrays,
		0, 0, 0, 0, 0};
	struct timeval start, end;
	gettimeofday(&start,NULL);
	swInt deIdx;
	for(deIdx=0;deIdx<dimension;deIdx++)
	{
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
			MLBFunParas.k2    = deIdx*edgeNumber;
			MLBFunParas.k3    = deIdx*vertexNumber;
			MLBParas->operatorFunPointer_host(&MLBFunParas);
//		for(k=startIdx;k<endIdx;k++)
//		{
//			b[neighbor[k]] += lower[k] * x[owner[k]];
//			b[owner[k]]    += upper[k] * x[neighbor[k]];
//		}
		}
	}

	for(spIndex=0;spIndex<mshBlockNum;spIndex++)
	{
		initTable(spIndex);
		e2VParas.spIndex = spIndex;
		athread_spawn(edge2VertexIter_slave, &e2VParas);
//		idxMP = spIndex+1;
//		for(deIdx=0;deIdx<dimension;deIdx++)
//		{
//			for(j=0;j<BLOCKNUM64K;j++)
//			{
//				row      = idxMP*BLOCKNUM64K+j;
//				col      = (idxMP+1)*BLOCKNUM64K;
////			col      = row;
//				if(col>=cpeBlockNum) break;
//				blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
//				startIdx = blockStarts[4*blockIdx+2];
//				col      = cpeBlockNum-1;
//				blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
//				endIdx   = blockStarts[4*blockIdx+3];
//				MLBFunParas.count = endIdx-startIdx;
//		    	MLBFunParas.k1    = startIdx;
//				MLBFunParas.k2    = deIdx*edgeNumber;
//				MLBFunParas.k3    = deIdx*vertexNumber;
//			    MLBParas->operatorFunPointer_host(&MLBFunParas);
////			for(k=startIdx;k<endIdx;k++)
////			{
////				b[neighbor[k]] +=lower[k] * x[owner[k]];
////				b[owner[k]]    +=upper[k] * x[neighbor[k]];
////			}
//			}
//		}
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

void initOwnNeiSendList(MLBParameters *MLBParas, swInt *ownNeiSendList){
	swInt  mshBlockNum = MLBParas->mshBlockNum;
	swInt  cpeBlockNum = MLBParas->cpeBlockNum;
	swInt* blockStarts = MLBParas->blockStarts;
	int num = mshBlockNum*BLOCKNUM64K*BLOCKNUM64K;
	int startIdx,endIdx,startListIdx,i,j,k,row,col,blockIdx,upperIdx,lowerIdx;
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

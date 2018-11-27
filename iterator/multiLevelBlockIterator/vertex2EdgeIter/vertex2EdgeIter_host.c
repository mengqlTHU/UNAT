#include "vertex2EdgeIter_host.h"
#include <stdio.h>

swInt spIndex, cpeBlockNum, mshBlockNum, mtxBlockNum;
swInt *blockStarts, *cellStarts, *ownNeiSendList;
swInt *firstEdgeVertices, *vertexNeighbor;
swFloat *data, *x, *b, *diag;
void (*operatorFunPointer_h)(MLBFunParameters *MLBFunParas);
void (*operatorFunPointer_s)(MLBFunParameters *MLBFunParas);

void vertex2EdgeIteration_host(Arrays* neighbourData, Arrays* vertexData,
			void(*operatorFunPointer_host) (MLBFunParameters *MLBFunParas),
			void(*operatorFunPointer_slave) (MLBFunParameters *MLBFunParas),
			MLBParameters* MLBParas)
{
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

	int i,j,k,row,col,startIdx,endIdx,blockIdx;

	for(i=0;i<cpeBlockNum;i++)
	{
		for(j=0;j<cpeBlockNum;j++)
		{
			blockIdx = i*cpeBlockNum+j;

		}
	}
	for(i=0;i<vertexData->num;i++)
	{
		b[i] = diag[i]*x[i];
	}
	for(i=0;i<neighbourData->num*2;i++)
	{
		b[firstEdgeVertices[i]] += data[i]*x[vertexNeighbor[i]];
	}
	printf("%f\n",b[30527]);

	LOG("edge2VertexIteration_host");

}


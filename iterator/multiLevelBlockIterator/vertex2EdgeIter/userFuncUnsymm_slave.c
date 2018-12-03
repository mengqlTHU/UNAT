// flag = 0 : 对角元素
// flag = 1 : 非对角元素
// integrate: x需要指定为NULL

#include "userFuncUnsymm_slave.h"
#include <stdlib.h>
#include "slave.h"

void spMVUnsymm(MLBFunParameters *MLBFunParas)
{
	swFloat* b               = MLBFunParas->vertexData->A1Ptr;
	swFloat* sx              = MLBFunParas->vertexData->A2Ptr;
	swFloat* rx              = MLBFunParas->vertexData->A3Ptr;
	swFloat* A1Ptr           = MLBFunParas->edgeData->A1Ptr;
	swInt* firstEdgeVertices = MLBFunParas->tArrays->sOwner;
	swInt* vertexNeighbor    = MLBFunParas->tArrays->sNeighbor;
	swInt  count         = MLBFunParas->count;
	swInt  k1            = MLBFunParas->k1;
	swInt  flag          = MLBFunParas->flag;
	int i;
	if(flag==0)
	{
		for(i=0;i<count;i++)
		{
			b[i] += sx[i]*rx[i];
		}
	} else if(flag==1)
	{
		for(i=0;i<count;i++)
		{
			b[firstEdgeVertices[i]-k1] += A1Ptr[i]*rx[i];
		}
	}
}

void integrateUnsymm(MLBFunParameters *MLBFunParas)
{
	swFloat* b               = MLBFunParas->vertexData->A1Ptr;
	swFloat* A1Ptr           = MLBFunParas->edgeData->A1Ptr;
	swInt* firstEdgeVertices = MLBFunParas->tArrays->sOwner;
	swInt  count         = MLBFunParas->count;
	swInt  k1            = MLBFunParas->k1;
	int i;
	for(i=0;i<count;i++)
	{
		b[firstEdgeVertices[i]-k1] += A1Ptr[i];
	}
}

// flag = 0 : 下三角块
// flag = 1 : 上三角块
// flag = 2 : 对角块上三角元素
// flag = 3 : 对角块下三角元素
// flag = 4 : 对角块对角元素
// flag = 5 : 主核稀疏块

#include "userFuncUnsymm_host.h"
#include <stdlib.h>
#include "slave.h"

void spMVUnsymm(MLBFunParameters *MLBFunParas)
{
	swFloat* b               = MLBFunParas->vertexData->A1Ptr;
	swFloat* x               = MLBFunParas->vertexData->A2Ptr;
	swFloat* data            = MLBFunParas->edgeData->A1Ptr;
	swInt* firstEdgeVertices = MLBFunParas->tArrays->sOwner;
	swInt* vertexNeighbor    = MLBFunParas->tArrays->sNeighbor;
	// count: endIdx-startIdx
	swInt  count             = MLBFunParas->count;
	//k1: startIdx
	swInt  k1                = MLBFunParas->k1;
	int i;
	for(i=k1;i<count+k1;i++)
	{
		b[firstEdgeVertices[i]] += data[i]*x[vertexNeighbor[i]];
	}
}

void integrateUnsymm(MLBFunParameters *MLBFunParas)
{
	swFloat* b               = MLBFunParas->vertexData->A1Ptr;
	swFloat* data            = MLBFunParas->edgeData->A1Ptr;
	swInt* firstEdgeVertices = MLBFunParas->tArrays->sOwner;
	swInt  count             = MLBFunParas->count;
	swInt  k1                = MLBFunParas->k1;
	int i;
	for(i=k1;i<count+k1;i++)
	{
		b[firstEdgeVertices[i]] += data[i];
	}

}

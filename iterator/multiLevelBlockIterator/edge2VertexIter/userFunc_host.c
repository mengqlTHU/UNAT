// flag = 0 : 下三角块
// flag = 1 : 上三角块
// flag = 2 : 对角块上三角元素
// flag = 3 : 对角块下三角元素
// flag = 4 : 对角块对角元素
// flag = 5 : 主核稀疏块

#include "userFunc_host.h"
#include <stdlib.h>
#include "slave.h"

void spMV(MLBFunParameters *MLBFunParas)
{
	swFloat* b           = MLBFunParas->vertexData->A1Ptr;
	swFloat* sx          = MLBFunParas->vertexData->A2Ptr;
	swFloat* A1Ptr       = MLBFunParas->edgeData->A1Ptr;
	swFloat* A2Ptr       = MLBFunParas->edgeData->A2Ptr;
	swInt* sOwner        = MLBFunParas->tArrays->sOwner;
	swInt* sNeighbor     = MLBFunParas->tArrays->sNeighbor;
	swInt  count         = MLBFunParas->count;
	swInt  k1            = MLBFunParas->k1;
	swInt  edgeNumber    = MLBFunParas->k2;
	swInt  vertexNumber  = MLBFunParas->k3;
	int i;
	for(i=k1;i<count+k1;i++)
	{
		b[sNeighbor[i]+vertexNumber]
			+= A1Ptr[i+edgeNumber]*sx[sOwner[i]+vertexNumber];
		b[sOwner[i]+vertexNumber]
			+= A2Ptr[i+edgeNumber]*sx[sNeighbor[i]+vertexNumber];
	}
}

void integrate(MLBFunParameters *MLBFunParas)
{
	swFloat* b           = MLBFunParas->vertexData->A1Ptr;
	swFloat* A1Ptr       = MLBFunParas->edgeData->A1Ptr;
	swFloat* A2Ptr       = MLBFunParas->edgeData->A2Ptr;
	swInt* sOwner        = MLBFunParas->tArrays->sOwner;
	swInt* sNeighbor     = MLBFunParas->tArrays->sNeighbor;
	swInt  count         = MLBFunParas->count;
	swInt  k1            = MLBFunParas->k1;
	int i;
	for(i=k1;i<count+k1;i++)
	{
		b[sNeighbor[i]] += A1Ptr[i];
		b[sOwner[i]]    += A2Ptr[i];
	}

}

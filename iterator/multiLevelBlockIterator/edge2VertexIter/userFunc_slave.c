// flag = 0 : 下三角块
// flag = 1 : 上三角块
// flag = 2 : 对角块上三角元素
// flag = 3 : 对角块下三角元素
// flag = 4 : 对角块对角元素
// flag = 5 : 主核稀疏块

#include "userFunc_slave.h"
#include <stdlib.h>
#include "slave.h"

void spMV(MLBFunParameters *MLBFunParas)
{
	swFloat* b           = MLBFunParas->vertexData->A1Ptr;
	swFloat* sx          = MLBFunParas->vertexData->A2Ptr;
	swFloat* rx          = MLBFunParas->vertexData->A3Ptr;
	swFloat* A1Ptr       = MLBFunParas->edgeData->A1Ptr;
	swFloat* A2Ptr       = MLBFunParas->edgeData->A2Ptr;
	swFloat* A3Ptr       = MLBFunParas->edgeData->A3Ptr;
	swInt* sOwner        = MLBFunParas->tArrays->sOwner;
	swInt* rOwner        = MLBFunParas->tArrays->rOwner;
	swInt* diagOwner     = MLBFunParas->tArrays->diagOwner;
	swInt* sNeighbor     = MLBFunParas->tArrays->sNeighbor;
	swInt* rNeighbor     = MLBFunParas->tArrays->rNeighbor;
	swInt* diagNeighbor  = MLBFunParas->tArrays->diagNeighbor;
	swInt  count         = MLBFunParas->count;
	swInt  k1            = MLBFunParas->k1;
	swInt  k2            = MLBFunParas->k2;
	swInt  flag          = MLBFunParas->flag;
	int i;
	if(flag == 0)
	{
		for(i=0;i<count;i++)
		{
			b[rNeighbor[i]-k1] += A2Ptr[i];
		}
	} else if(flag == 1)
	{
		for(i=0;i<count;i++)
		{
			b[sNeighbor[i]-k1] += A1Ptr[i]*rx[i];
		}
	} else if(flag == 2)
	{
		for(i=0;i<count;i++)
		{
			b[diagOwner[i]-k1] += A3Ptr[i]*sx[diagNeighbor[i]-k1];
		}
	} else if(flag == 3)
	{
		for(i=0;i<count;i++)
		{
			b[diagNeighbor[i]-k1] += A3Ptr[i]*sx[diagOwner[i]-k1];
		}
	} else if(flag == 4)
	{
		for(i=0;i<count;i++)
		{
			b[i] += rx[i]*sx[i];
		}
	} else if(flag == 5)
	{
		LOG("test");
//		for(i=k1;i<count+k1;i++)
//		{
//			b[sNeighbor[i]] += A1Ptr[i]*sx[sOwner[i]];
//			b[sOwner[i]]    += A2Ptr[i]*sx[sNeighbor[i]];
//		}

	} else
	{
		LOG("The flag is not defined");
	}
}

void integrate(MLBFunParameters *MLBFunParas)
{
	swFloat* b           = MLBFunParas->vertexData->A1Ptr;
	swFloat* sx          = MLBFunParas->vertexData->A2Ptr;
	swFloat* rx          = MLBFunParas->vertexData->A3Ptr;
	swFloat* A1Ptr       = MLBFunParas->edgeData->A1Ptr;
	swFloat* A2Ptr       = MLBFunParas->edgeData->A2Ptr;
	swFloat* A3Ptr       = MLBFunParas->edgeData->A3Ptr;
	swInt* sOwner        = MLBFunParas->tArrays->sOwner;
	swInt* rOwner        = MLBFunParas->tArrays->rOwner;
	swInt* diagOwner     = MLBFunParas->tArrays->diagOwner;
	swInt* sNeighbor     = MLBFunParas->tArrays->sNeighbor;
	swInt* rNeighbor     = MLBFunParas->tArrays->rNeighbor;
	swInt* diagNeighbor  = MLBFunParas->tArrays->diagNeighbor;
	swInt  count         = MLBFunParas->count;
	swInt  k1            = MLBFunParas->k1;
	swInt  k2            = MLBFunParas->k2;
	swInt  flag          = MLBFunParas->flag;
	int i;
	if(flag == 0)
	{
		for(i=0;i<count;i++)
		{
			b[rNeighbor[i]-k1] += A2Ptr[i];
		}
	} else if(flag == 1)
	{
		for(i=0;i<count;i++)
		{
			b[sNeighbor[i]-k1] += A1Ptr[i];
		}
	} else if(flag == 2)
	{
		for(i=0;i<count;i++)
		{
			b[diagOwner[i]-k1] += A3Ptr[i];
		}
	} else if(flag == 3)
	{
		for(i=0;i<count;i++)
		{
			b[diagNeighbor[i]-k1] += A3Ptr[i];
		}
	} else if(flag == 4)
	{
		for(i=0;i<count;i++)
		{
			b[i] += rx[i]*sx[i];
		}
	} else
	{
		LOG("The flag is not defined");
	}
}

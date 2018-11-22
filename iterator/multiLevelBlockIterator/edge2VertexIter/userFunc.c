#include "userFunc.h"
#include <stdlib.h>
#include "slave.h"

void spMV(MLBFunParameters *MLBFunParas)
{
	swFloat* b     = MLBFunParas->vertexData->A1Ptr;
	swFloat* x     = MLBFunParas->vertexData->A2Ptr;
	swFloat* diag  = MLBFunParas->vertexData->A3Ptr;
	swFloat* lower = MLBFunParas->edgeData->A1Ptr;
	swFloat* upper = MLBFunParas->edgeData->A2Ptr;
	swInt* owner   = MLBFunParas->owner;
	swInt* neighbor= MLBFunParas->neighbor;
	swInt  k1      = MLBFunParas->k1;
	swInt  k2      = MLBFunParas->k2;
	int i;
	for(i=0;i<MLBFunParas->count;i++)
	{
		b[neighbor[i]-k1] += lower[i]*x[i];
	}
}

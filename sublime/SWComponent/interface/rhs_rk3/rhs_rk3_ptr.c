#include "swMacro.h"
#include "iterator.h"
#include "math.h"

define_e2v_FunPtr(rhs_rk3_upt_rhs)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_rhs  = accessArray(data->frontEdgeData, 0);
	swFloat* fv_duc  = accessArray(data->frontEdgeData, 1);
	swFloat* fv_rv0  = accessArray(data->frontEdgeData, 2);

	swFloat* para = accessArray(paraData, 0);

	swInt edgeNum = getArraySize(data->frontEdgeData);

	for(iedge=0;iedge<edgeNum;iedge++)
	{
		if(para[0]==1)
		{
			fv_rv0[iedge*5+0] = fv_duc[iedge*5+0];
			fv_rv0[iedge*5+1] = fv_duc[iedge*5+1];
			fv_rv0[iedge*5+2] = fv_duc[iedge*5+2];
			fv_rv0[iedge*5+3] = fv_duc[iedge*5+3];
			fv_rv0[iedge*5+4] = fv_duc[iedge*5+4];
		} else
		{
			fv_duc[iedge*5+0] = 0.5*(fv_rv0[iedge*5+0]+fv_duc[iedge*5+0]);
			fv_duc[iedge*5+1] = 0.5*(fv_rv0[iedge*5+1]+fv_duc[iedge*5+1]);
			fv_duc[iedge*5+2] = 0.5*(fv_rv0[iedge*5+2]+fv_duc[iedge*5+2]);
			fv_duc[iedge*5+3] = 0.5*(fv_rv0[iedge*5+3]+fv_duc[iedge*5+3]);
			fv_duc[iedge*5+4] = 0.5*(fv_rv0[iedge*5+4]+fv_duc[iedge*5+4]);
		}
		fv_rhs[iedge*5+0] -= fv_duc[iedge*5+0];
		fv_rhs[iedge*5+1] -= fv_duc[iedge*5+1];
		fv_rhs[iedge*5+2] -= fv_duc[iedge*5+2];
		fv_rhs[iedge*5+3] -= fv_duc[iedge*5+3];
		fv_rhs[iedge*5+4] -= fv_duc[iedge*5+4];
	}
}

define_e2v_FunPtr(rhs_rk3_zero)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_duc  = accessArray(data->frontEdgeData, 0);
	swFloat* fv_rhs  = accessArray(data->frontEdgeData, 1);
	swInt edgeNum = getArraySize(data->frontEdgeData);
	for(iedge=0;iedge<edgeNum*5;iedge++)
	{
		fv_duc[iedge] = 0.0;
		fv_rhs[iedge] = 0.0;
	}
}

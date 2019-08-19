#include "swMacro.h"
#include "iterator.h"
#include "math.h"

define_e2v_FunPtr(exp_rk43_copy)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_uc0  = accessArray(data->frontEdgeData, 0);
	swInt dims_uc0   = getArrayDims(data->frontEdgeData, 0);
	swFloat* fv_uc   = accessArray(data->frontEdgeData, 1);
	swInt edgeNum = getArraySize(data->frontEdgeData);

	for(iedge=0;iedge<edgeNum*dims_uc0;iedge++)
	{
		fv_uc0[iedge] = fv_uc[iedge];
	}
}

define_e2v_FunPtr(exp_rk43_1st)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_uc    = accessArray(data->frontEdgeData, 0);
	swInt dims_uc     = getArrayDims(data->frontEdgeData, 0);
	swFloat* fv_uc0   = accessArray(data->frontEdgeData, 1);
	swInt dims_uc0    = getArrayDims(data->frontEdgeData, 1);
	swFloat* fv_rhs   = accessArray(data->frontEdgeData, 2);
	swInt dims_rhs    = getArrayDims(data->frontEdgeData, 2);
	swFloat* fv_LHS_s = accessArray(data->frontEdgeData, 3);
	swFloat* is_int   = accessArray(data->frontEdgeData, 4);

	swInt edgeNum = getArraySize(data->frontEdgeData);

	for(iedge=0;iedge<edgeNum;iedge++)
	{
		if(is_int[iedge]==1)
		{
			fv_uc[iedge*dims_uc+0] = fv_uc0[iedge*dims_uc0+0]
				- 0.5*fv_rhs[iedge*dims_rhs+0]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+1] = fv_uc0[iedge*dims_uc0+1]
				- 0.5*fv_rhs[iedge*dims_rhs+1]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+2] = fv_uc0[iedge*dims_uc0+2]
				- 0.5*fv_rhs[iedge*dims_rhs+2]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+3] = fv_uc0[iedge*dims_uc0+3]
				- 0.5*fv_rhs[iedge*dims_rhs+3]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+4] = fv_uc0[iedge*dims_uc0+4]
				- 0.5*fv_rhs[iedge*dims_rhs+4]*fv_LHS_s[iedge];
		}
	}
}

define_e2v_FunPtr(exp_rk43_2nd)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_uc    = accessArray(data->frontEdgeData, 0);
	swInt dims_uc     = getArrayDims(data->frontEdgeData, 0);
	swFloat* fv_rhs   = accessArray(data->frontEdgeData, 1);
	swInt dims_rhs    = getArrayDims(data->frontEdgeData, 1);
	swFloat* fv_LHS_s = accessArray(data->frontEdgeData, 2);
	swFloat* is_int   = accessArray(data->frontEdgeData, 3);

	swInt edgeNum = getArraySize(data->frontEdgeData);

	for(iedge=0;iedge<edgeNum;iedge++)
	{
		if(is_int[iedge]==1)
		{
			fv_uc[iedge*dims_uc+0] -=
				0.5*fv_rhs[iedge*dims_rhs+0]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+1] -=
				0.5*fv_rhs[iedge*dims_rhs+1]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+2] -=
				0.5*fv_rhs[iedge*dims_rhs+2]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+3] -= 
				0.5*fv_rhs[iedge*dims_rhs+3]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+4] -=
				0.5*fv_rhs[iedge*dims_rhs+4]*fv_LHS_s[iedge];
		}
	}
}

define_e2v_FunPtr(exp_rk43_3rd)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_uc    = accessArray(data->frontEdgeData, 0);
	swInt dims_uc     = getArrayDims(data->frontEdgeData, 0);
	swFloat* fv_uc0   = accessArray(data->frontEdgeData, 1);
	swInt dims_uc0    = getArrayDims(data->frontEdgeData, 1);
	swFloat* fv_rhs   = accessArray(data->frontEdgeData, 2);
	swInt dims_rhs    = getArrayDims(data->frontEdgeData, 2);
	swFloat* fv_LHS_s = accessArray(data->frontEdgeData, 3);
	swFloat* is_int   = accessArray(data->frontEdgeData, 4);

	swFloat* para = accessArray(paraData, 0);

	swInt edgeNum = getArraySize(data->frontEdgeData);

	for(iedge=0;iedge<edgeNum;iedge++)
	{
		if(is_int[iedge]==1)
		{
			fv_uc[iedge*dims_uc+0] = para[0]*fv_uc0[iedge*dims_uc0+0]
				+ para[1]*fv_uc[iedge*dims_uc+0]
				- para[2]*fv_rhs[iedge*dims_rhs+0]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+1] = para[0]*fv_uc0[iedge*dims_uc0+1]
				+ para[1]*fv_uc[iedge*dims_uc+1]
				- para[2]*fv_rhs[iedge*dims_rhs+1]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+2] = para[0]*fv_uc0[iedge*dims_uc0+2]
				+ para[1]*fv_uc[iedge*dims_uc+2]
				- para[2]*fv_rhs[iedge*dims_rhs+2]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+3] = para[0]*fv_uc0[iedge*dims_uc0+3]
				+ para[1]*fv_uc[iedge*dims_uc+3]
				- para[2]*fv_rhs[iedge*dims_rhs+3]*fv_LHS_s[iedge];
			fv_uc[iedge*dims_uc+4] = para[0]*fv_uc0[iedge*dims_uc0+4]
				+ para[1]*fv_uc[iedge*dims_uc+4]
				- para[2]*fv_rhs[iedge*dims_rhs+4]*fv_LHS_s[iedge];
		}
	}
}

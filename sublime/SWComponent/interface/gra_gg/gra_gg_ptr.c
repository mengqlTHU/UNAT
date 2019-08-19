#include "swMacro.h"
#include "iterator.h"
#include "math.h"

//define 2 function pointers
define_e2v_FunPtr(gra_gg_div)
{
	//selfConn computation
	int iedge,j,idx;
	swFloat* fv_gra  = accessArray(data->frontEdgeData, 0);
	swInt dims_gra   = getArrayDims(data->frontEdgeData, 0);
	swFloat* sec_vol = accessArray(data->frontEdgeData, 1);
	swInt edgeNum = getArraySize(data->frontEdgeData);
	for(iedge=0;iedge<edgeNum;iedge++)
	{
		idx = iedge*dims_gra;
		for(j=0;j<18;j++)
		{
			fv_gra[idx+j] /= sec_vol[iedge];
		}
	}
}

define_e2v_FunPtr(gra_gg_zero)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_gra  = accessArray(data->frontEdgeData, 0);
	swInt dims_gra   = getArrayDims(data->frontEdgeData, 0);
	swFloat* sec_vol = accessArray(data->frontEdgeData, 1);
	swInt edgeNum = getArraySize(data->frontEdgeData);
	for(iedge=0;iedge<edgeNum*dims_gra;iedge++)
	{
		fv_gra[iedge] = 0.0;
	}
}

#include "swMacro.h"
#include "iterator.h"
#include "math.h"

//define 2 function pointers
define_e2v_FunPtr(uc2u)
{
	//selfConn computation
	int iedge,j,idx;
// printf("hello world\n");
	swFloat* fv_u   = accessArray(data->frontEdgeData, 0);
	swInt dims_u    = getArrayDims(data->frontEdgeData, 0);
	swFloat* fv_uc  = accessArray(data->frontEdgeData, 1);
	swInt dims_uc   = getArrayDims(data->frontEdgeData, 1);
	swFloat* fv_t   = accessArray(data->frontEdgeData, 2);
	swFloat* is_int = accessArray(data->frontEdgeData, 3);
	swFloat* fv_mu  = accessArray(data->frontEdgeData, 4);
	swInt dims_mu   = getArrayDims(data->frontEdgeData, 4);
	swFloat* para   = accessArray(paraData, 0);

	swInt edgeNum = getArraySize(data->frontEdgeData);
	swFloat ke;
	for(iedge=0;iedge<edgeNum;iedge++)
	{
		idx = iedge*dims_u;
		if(is_int[iedge]==1)
		{
			fv_u[idx+0] = fv_uc[idx+0];
			fv_u[idx+1] = fv_uc[idx+1]/fv_uc[idx+0];
			fv_u[idx+2] = fv_uc[idx+2]/fv_uc[idx+0];
			fv_u[idx+3] = fv_uc[idx+3]/fv_uc[idx+0];
			ke = 0.5*(fv_u[idx+1]*fv_u[idx+1]
				+ fv_u[idx+2]*fv_u[idx+2]
				+ fv_u[idx+3]*fv_u[idx+3]);
			fv_u[idx+4] 
				= fabs(fv_uc[idx+4]-fv_uc[idx+0]*ke)*para[1];
			fv_t[iedge] = fv_u[idx+4]/(para[0]*fv_uc[idx+0]);
			fv_mu[iedge*dims_mu+0] = 1.461e-6*sqrt(fv_t[iedge]*fv_t[iedge]*fv_t[iedge])/
				(fv_t[iedge]+para[2]);
		}
	}
}

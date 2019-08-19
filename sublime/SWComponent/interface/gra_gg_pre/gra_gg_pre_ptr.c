#include "swMacro.h"
#include "iterator.h"
#include "math.h"

// #ifdef __cplusplus
// extern "C"
// {
// #endif

//define 2 function pointers
define_e2v_FunPtr(gra_gg_pre)
{
	//selfConn computation
	int i,j;
	swFloat* fv_gra    = accessArray(data->vertexData, 0);
	swInt dims_gra     = getArrayDims(data->vertexData, 0);
	swFloat* fv_u      = accessArray(data->vertexData, 1);
	swInt dims_u       = getArrayDims(data->vertexData, 1);
	swFloat* sec_vol   = accessArray(data->vertexData, 2);
	swFloat* fv_t      = accessArray(data->vertexData, 3);
	swFloat* is_int    = accessArray(data->vertexData, 4);
	
	//frontEdge computation
	swFloat* n_vg	      = accessArray(data->frontEdgeData, 0);
	swInt dims_vg         = getArrayDims(data->frontEdgeData, 0);
	swFloat* transform    = accessArray(data->frontEdgeData, 1);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt iedge,idxL,idxR;
	swFloat vL,vR,um[6],n[3];
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		if(transform[iedge] == 0)
		{
			idxL = startVertices[iedge];
			idxR = endVertices[iedge];
		} else
		{
			idxL = endVertices[iedge];
			idxR = startVertices[iedge];
		}
		n[0] = n_vg[iedge*dims_vg+0]*n_vg[iedge*dims_vg+3];
		n[1] = n_vg[iedge*dims_vg+1]*n_vg[iedge*dims_vg+3];
		n[2] = n_vg[iedge*dims_vg+2]*n_vg[iedge*dims_vg+3];

		vL = 1.0/sec_vol[idxL];
		vR = 1.0/sec_vol[idxR];
		um[0] = (vL*fv_u[idxL*dims_u+0]+vR*fv_u[idxR*dims_u+0])/(vL+vR);
		um[1] = (vL*fv_u[idxL*dims_u+1]+vR*fv_u[idxR*dims_u+1])/(vL+vR);
		um[2] = (vL*fv_u[idxL*dims_u+2]+vR*fv_u[idxR*dims_u+2])/(vL+vR);
		um[3] = (vL*fv_u[idxL*dims_u+3]+vR*fv_u[idxR*dims_u+3])/(vL+vR);
		um[4] = (vL*fv_u[idxL*dims_u+4]+vR*fv_u[idxR*dims_u+4])/(vL+vR);
		um[5] = (vL*fv_t[idxL         ]+vR*fv_t[idxR         ])/(vL+vR);
		// fv_gra[idxL*dims_gra+0] += um[5];
		// fv_gra[idxR*dims_gra+0] += um[0];
		for(i=0;i<6;i++)
		{
			fv_gra[idxL*dims_gra+3*i  ] += um[i]*n[0];
			fv_gra[idxL*dims_gra+3*i+1] += um[i]*n[1];
			fv_gra[idxL*dims_gra+3*i+2] += um[i]*n[2];
		}
		if(is_int[idxR]==1)
		{
			for(i=0;i<6;i++)
			{
				fv_gra[idxR*dims_gra+3*i  ] -= um[i]*n[0];
				fv_gra[idxR*dims_gra+3*i+1] -= um[i]*n[1];
				fv_gra[idxR*dims_gra+3*i+2] -= um[i]*n[2];
			}			
		}
		// printf("%d\n", dims_gra);
	}
}

// #ifdef __cplusplus
// }
// #endif

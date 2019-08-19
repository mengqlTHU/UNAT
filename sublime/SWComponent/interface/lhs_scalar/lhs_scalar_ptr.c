#include "swMacro.h"
#include "iterator.h"
#include "math.h"


//define 2 function pointers
define_e2v_FunPtr(lhs_scalar)
{
	//selfConn computation
	int i,j;
	swFloat* fv_LHS_s = accessArray(data->vertexData, 0);
	swInt dims_lhs_s  = getArrayDims(data->vertexData, 0);
	swFloat* fv_u     = accessArray(data->vertexData, 1);
	swInt dims_u      = getArrayDims(data->vertexData, 1);
	swFloat* is_int   = accessArray(data->vertexData, 2);
	swFloat* fv_mu    = accessArray(data->vertexData, 3);
	swInt dims_mu     = getArrayDims(data->vertexData, 3);
	
	//frontEdge computation
	swFloat* n_vg	      = accessArray(data->frontEdgeData, 0);
	swInt dims_vg         = getArrayDims(data->frontEdgeData, 0);
	swFloat* transform    = accessArray(data->frontEdgeData, 1);

	swFloat* para = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt iedge,idxL,idxR;
	swFloat mu,mut,um[6],u[5],a,s;
	swFloat *uL,*uR,*n;
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
		uL   = &fv_u[idxL*dims_u+0];
		uR   = &fv_u[idxR*dims_u+0];
		u[0] = 0.5*(uL[0]+uR[0]);
		u[1] = 0.5*(uL[1]+uR[1]);
		u[2] = 0.5*(uL[2]+uR[2]);
		u[3] = 0.5*(uL[3]+uR[3]);
		u[4] = 0.5*(uL[4]+uR[4]);
		if(para[0]==1) mu = 0.5*(fv_mu[idxL*dims_mu+0]+fv_mu[idxR*dims_mu+0]);
		if(para[1]==1) mut = 0.5*(fv_mu[idxL*dims_mu+1]+fv_mu[idxR*dims_mu+1]);
		a = sqrt(para[3]*u[4]/u[0]);

		n = &n_vg[iedge*dims_vg+0];

		u[1] -= n_vg[iedge*dims_vg+5];
		u[2] -= n_vg[iedge*dims_vg+6];
		u[3] -= n_vg[iedge*dims_vg+7];

		s =  0.5*n[3]*(fabs(u[1]*n[0]+u[2]*n[1]+u[3]*n[2])+a);

		if(para[2]==-1)
		{
            if(para[0]==1)  s += para[3]*mu *n[3]*n[4]/u[0];
            if(para[1]==1)  s += para[3]*mut*n[3]*n[4]/u[0];
		}

		fv_LHS_s[idxL] += s;
		if(is_int[idxR]==1) fv_LHS_s[idxR] += s;
	}
}

define_e2v_FunPtr(lhs_scalar_zero)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_LHS_s  = accessArray(data->frontEdgeData, 0);
	swInt edgeNum = getArraySize(data->frontEdgeData);
	for(iedge=0;iedge<edgeNum;iedge++)
	{
		fv_LHS_s[iedge] = 0.0;
	}
}

define_e2v_FunPtr(lhs_scalar_cal_dt)
{
	//selfConn computation
	int iedge,j;
	swFloat* fv_LHS_s  = accessArray(data->frontEdgeData, 0);
	swFloat* fv_uc     = accessArray(data->frontEdgeData, 1);
	swInt dims_uc      = getArrayDims(data->frontEdgeData, 1);
	swFloat* fv_rhs    = accessArray(data->frontEdgeData, 2);
	swInt dims_rhs     = getArrayDims(data->frontEdgeData, 2);
	swFloat* is_int    = accessArray(data->frontEdgeData, 3);

	swFloat* para = accessArray(paraData, 0);

	swInt edgeNum = getArraySize(data->frontEdgeData);
	swFloat s,CFLu;
	swFloat *uc,*R;
	for(iedge=0;iedge<edgeNum;iedge++)
	{
		if(is_int[iedge]==1)
		{
			s  = fv_LHS_s[iedge];
			uc = &fv_uc[iedge*dims_uc+0];
			R  = &fv_rhs[iedge*dims_rhs+0];

			CFLu = para[0];
			CFLu = MIN(CFLu, para[1]*uc[0]*s/(fabs(R[0])+1.0e-9));
			CFLu = MIN(CFLu, para[1]*uc[4]*s/(fabs(R[4])+1.0e-9));
			CFLu = MAX(CFLu, 1.0e-1*para[0]);

			fv_LHS_s[iedge] = CFLu/s;
		}
	}
}

#include "swMacro.h"
#include "iterator.h"
#include "math.h"

#define get_vis_fac(M, e, u, mu, g, k, rhsD) \
{ \
	swFloat mul, mut, kq, dU3, S[9],Tlam[9], Ttur[9], f[4], turk; \
	mul = mu[0]; \
	dU3 = (g[0]+g[4]+g[8])/3.0; \
	S[0] = g[0]-dU3; \
	S[1] = 0.5*(g[1]+g[3]); \
	S[2] = 0.5*(g[2]+g[6]); \
	S[3] = S[1]; \
	S[4] = g[4]-dU3; \
	S[5] = 0.5*(g[5]+g[7]); \
	S[6] = S[2]; \
	S[7] = S[5]; \
	S[8] = g[8]-dU3; \
	swInt i_macro; \
	for(i_macro=0;i_macro<9;i_macro++) \
	{ \
		Tlam[i_macro] = 2.0*mul*S[i_macro]; \
	} \
	kq   = para[2]*mul/para[3]; \
	\
	if(para[1]==1) \
	{ \
		mut = mu[1]; \
		for(i_macro=0;i_macro<9;i_macro++) \
		{ \
			Ttur[i_macro] = 2.0*mut*S[i_macro]; \
		} \
		if(para[5]==1) printf("Error: is_KO has not been tested.\n"); \
		kq = kq+para[2]*mut/para[4]; \
		for(i_macro=0;i_macro<9;i_macro++) \
		{ \
			Tlam[i_macro] += Ttur[i_macro]; \
		} \
	} \
	\
	f[0] = e[0]*Tlam[0]+e[1]*Tlam[3]+e[2]*Tlam[6]; \
	f[1] = e[0]*Tlam[1]+e[1]*Tlam[4]+e[2]*Tlam[7]; \
	f[2] = e[0]*Tlam[2]+e[1]*Tlam[5]+e[2]*Tlam[8]; \
	f[3] = u[1]*f[0]+u[2]*f[1]+u[3]*f[2]+kq*(e[0]*g[9]+e[1]*g[10]+e[2]*g[11]); \
	rhsD[0] = 0.0; \
	rhsD[1] = f[0]*e[3]; \
	rhsD[2] = f[1]*e[3]; \
	rhsD[3] = f[2]*e[3]; \
	rhsD[4] = f[3]*e[3]; \
}

//define 2 function pointers
define_e2v_FunPtr(rhs_vis)
{
	//selfConn computation
	int i,j,k;
	swFloat* fv_rhs    = accessArray(data->vertexData, 0);
	swInt dims_rhs     = getArrayDims(data->vertexData, 0);
	swFloat* fv_u      = accessArray(data->vertexData, 1);
	swInt dims_u       = getArrayDims(data->vertexData, 1);
	swFloat* sec_cen   = accessArray(data->vertexData, 2);
	swInt dims_sec_cen = getArrayDims(data->vertexData, 2);
	swFloat* fv_gra    = accessArray(data->vertexData, 3);
	swInt dims_gra     = getArrayDims(data->vertexData, 3);
	swFloat* fv_duc    = accessArray(data->vertexData, 4);
	swInt dims_duc     = getArrayDims(data->vertexData, 4);
	swFloat* is_int    = accessArray(data->vertexData, 5);
	swFloat* fv_t      = accessArray(data->vertexData, 6);
	swFloat* fv_mu     = accessArray(data->vertexData, 7);
	swInt dims_mu      = getArrayDims(data->vertexData, 7);
	swFloat* fv_turb   = accessArray(data->vertexData, 8);
	swInt dims_turb    = getArrayDims(data->vertexData, 8);
	
	//frontEdge computation
	swFloat* n_vg	      = accessArray(data->frontEdgeData, 0);
	swInt dims_vg         = getArrayDims(data->frontEdgeData, 0);
	swFloat* transform    = accessArray(data->frontEdgeData, 1);
	swInt dims_trans      = getArrayDims(data->frontEdgeData, 1);

	swFloat* para = accessArray(paraData, 0);

	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt iedge,idxL,idxR;
	swFloat *fac_1d,*uLf,*uRf,g[12],u[5],d[3],tL,tR,du[4];
	swFloat nne[3],dr,mu[2],tke,rtmp,tmp,rhsD[5];
	for(iedge=0;iedge<edgeNumber;iedge++)
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
		fac_1d = &n_vg[iedge*dims_vg+0];

		uLf = &fv_u[idxL*dims_u+0];

		uRf = &fv_u[idxR*dims_u+0];

		u[0] = 0.5*(uLf[0]+uRf[0]);
		u[1] = 0.5*(uLf[1]+uRf[1]);
		u[2] = 0.5*(uLf[2]+uRf[2]);
		u[3] = 0.5*(uLf[3]+uRf[3]);
		u[4] = 0.5*(uLf[4]+uRf[4]);

		mu[0] = 0.5*(fv_mu[idxL*dims_mu+0]+fv_mu[idxR*dims_mu+0]);
		if(para[1]==1) mu[1] = 0.5*(fv_mu[idxL*dims_mu+1]+fv_mu[idxR*dims_mu+1]);
		if(para[5]==1) tke   = 0.5*(fv_turb[idxL*dims_turb+0]+fv_turb[idxR*dims_turb+0]);
		tL = fv_t[idxL];
		tR = fv_t[idxR];

		for(j=0;j<9;j++)
		{
			g[j] = 0.5*(fv_gra[idxL*dims_gra+j+3]+fv_gra[idxR*dims_gra+j+3]);
		}
		g[9]  = 0.5*(fv_gra[idxL*dims_gra+15]+fv_gra[idxR*dims_gra+15]);
		g[10] = 0.5*(fv_gra[idxL*dims_gra+16]+fv_gra[idxR*dims_gra+16]);
		g[11] = 0.5*(fv_gra[idxL*dims_gra+17]+fv_gra[idxR*dims_gra+17]);

		d[0]  = sec_cen[idxR*dims_sec_cen+0]-sec_cen[idxL*dims_sec_cen+0];
		d[1]  = sec_cen[idxR*dims_sec_cen+1]-sec_cen[idxL*dims_sec_cen+1];
		d[2]  = sec_cen[idxR*dims_sec_cen+2]-sec_cen[idxL*dims_sec_cen+2];
		rtmp  = 1.0/sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		d[0]  = d[0]*rtmp;
		d[1]  = d[1]*rtmp;
		d[2]  = d[2]*rtmp;
		du[0] = rtmp*(uRf[1]-uLf[1]);
		du[1] = rtmp*(uRf[2]-uLf[2]);
		du[2] = rtmp*(uRf[3]-uLf[3]);
		du[3] = rtmp*(tR-tL);

		tmp = fac_1d[0]*d[0]+fac_1d[1]*d[1]+fac_1d[2]*d[2];
		nne[0] = fac_1d[0]/tmp;
		nne[1] = fac_1d[1]/tmp;
		nne[2] = fac_1d[2]/tmp;

		for (j = 0; j < 4; ++j)
		{
			k       = 3*j;
			dr      = du[j]-g[k]*d[0]-g[k+1]*d[1]-g[k+2]*d[2];
			g[k]   += dr*nne[0];
			g[k+1] += dr*nne[1];
			g[k+2] += dr*nne[2];
		}

		get_vis_fac(1, fac_1d, u, mu, g, tke, rhsD);

		if(para[0]==1)
		{
			for(j=0;j<5;j++)
			{
				fv_rhs[idxL*dims_rhs+j] -= rhsD[j];
			}
			if (is_int[idxR]==1)
			{
				for(j=0;j<5;j++)
				{	
					fv_rhs[idxR*dims_rhs+j] += rhsD[j];
				}
			}
		} else
		{
			for(j=0;j<5;j++)
			{
				fv_duc[idxL*dims_rhs+j] += rhsD[j];
			}
			if (is_int[idxR]==1)
			{
				for(j=0;j<5;j++)
				{	
					fv_duc[idxR*dims_rhs+j] -= rhsD[j];
				}
			}			
		}
	}
}



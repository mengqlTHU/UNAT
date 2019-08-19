#include "swMacro.h"
#include "iterator.h"
#include "math.h"

#define reigl_new(u,vg,e,sgn_eig,sgn_mat,eps,RDL,spra,eig) \
{ \
	swFloat a,qn,vn,eiga,eigb,m[5],n[5],ke,H,a1,a2; \
    a       =  sqrt(para[2]*u[4]/u[0]); \
    qn      =  e[0]*u [1]+e[1]*u [2]+e[2]*u [3]; \
    vn      =  e[0]*vg[0]+e[1]*vg[1]+e[2]*vg[2]; \
    m[0]    =  qn-vn; \
    m[1]    =  m[0]+a; \
    m[2]    =  m[0]-a; \
    spra    =  fabs(m[0])+a; \
    eig[0]  =  0.5*sgn_mat*e[3]*(m[0]+sgn_eig*(fabs(m[0])+eps*spra)); \
    eig[1]  =  0.5*sgn_mat*e[3]*(m[1]+sgn_eig*(fabs(m[1])+eps*spra)); \
    eig[2]  =  0.5*sgn_mat*e[3]*(m[2]+sgn_eig*(fabs(m[2])+eps*spra)); \
    eiga    =  0.5*(eig[1]-eig[2]); \
    eigb    =  0.5*(eig[1]+eig[2])-eig[0]; \
    \
    m[0]    = -qn; \
    m[1]    =  e[0]; \
    m[2]    =  e[1]; \
    m[3]    =  e[2]; \
    m[4]    =  0.0; \
    ke      =  0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3]); \
    n[0]    =  para[3]*ke; \
    n[1]    = -para[3]*u[1]; \
    n[2]    = -para[3]*u[2]; \
    n[3]    = -para[3]*u[3]; \
    n[4]    =  para[3]; \
    H       =  para[2]*u[4]/(para[3]*u[0])+ke; \
    \
    a1          =  1.0/a; \
    a2          =  a1*a1; \
    RDL[0][0]   =  eiga*a1*      m[0]          +eigb                *a2*n[0]; \
    RDL[0][1]   =  eiga*a1*      m[1]          +eigb                *a2*n[1]; \
    RDL[0][2]   =  eiga*a1*      m[2]          +eigb                *a2*n[2]; \
    RDL[0][3]   =  eiga*a1*      m[3]          +eigb                *a2*n[3]; \
    RDL[0][4]   =  eiga*a1*      m[4]          +eigb                *a2*n[4]; \
    \
    RDL[1][0]   =  eiga*a1*(u[1]*m[0]+e[0]*n[0])+eigb*(e[0]*m[0]+u[1]*a2*n[0]); \
    RDL[1][1]   =  eiga*a1*(u[1]*m[1]+e[0]*n[1])+eigb*(e[0]*m[1]+u[1]*a2*n[1]); \
    RDL[1][2]   =  eiga*a1*(u[1]*m[2]+e[0]*n[2])+eigb*(e[0]*m[2]+u[1]*a2*n[2]); \
    RDL[1][3]   =  eiga*a1*(u[1]*m[3]+e[0]*n[3])+eigb*(e[0]*m[3]+u[1]*a2*n[3]); \
    RDL[1][4]   =  eiga*a1*(u[1]*m[4]+e[0]*n[4])+eigb*(e[0]*m[4]+u[1]*a2*n[4]); \
    \
    RDL[2][0]   =  eiga*a1*(u[2]*m[0]+e[1]*n[0])+eigb*(e[1]*m[0]+u[2]*a2*n[0]); \
    RDL[2][1]   =  eiga*a1*(u[2]*m[1]+e[1]*n[1])+eigb*(e[1]*m[1]+u[2]*a2*n[1]); \
    RDL[2][2]   =  eiga*a1*(u[2]*m[2]+e[1]*n[2])+eigb*(e[1]*m[2]+u[2]*a2*n[2]); \
    RDL[2][3]   =  eiga*a1*(u[2]*m[3]+e[1]*n[3])+eigb*(e[1]*m[3]+u[2]*a2*n[3]); \
    RDL[2][4]   =  eiga*a1*(u[2]*m[4]+e[1]*n[4])+eigb*(e[1]*m[4]+u[2]*a2*n[4]); \
    \
    RDL[3][0]   =  eiga*a1*(u[3]*m[0]+e[2]*n[0])+eigb*(e[2]*m[0]+u[3]*a2*n[0]); \
    RDL[3][1]   =  eiga*a1*(u[3]*m[1]+e[2]*n[1])+eigb*(e[2]*m[1]+u[3]*a2*n[1]); \
    RDL[3][2]   =  eiga*a1*(u[3]*m[2]+e[2]*n[2])+eigb*(e[2]*m[2]+u[3]*a2*n[2]); \
    RDL[3][3]   =  eiga*a1*(u[3]*m[3]+e[2]*n[3])+eigb*(e[2]*m[3]+u[3]*a2*n[3]); \
    RDL[3][4]   =  eiga*a1*(u[3]*m[4]+e[2]*n[4])+eigb*(e[2]*m[4]+u[3]*a2*n[4]); \
    \
    RDL[4][0]   =  eiga*a1*(H   *m[0]+qn  *n[0])+eigb*(qn  *m[0]+H   *a2*n[0]); \
    RDL[4][1]   =  eiga*a1*(H   *m[1]+qn  *n[1])+eigb*(qn  *m[1]+H   *a2*n[1]); \
    RDL[4][2]   =  eiga*a1*(H   *m[2]+qn  *n[2])+eigb*(qn  *m[2]+H   *a2*n[2]); \
    RDL[4][3]   =  eiga*a1*(H   *m[3]+qn  *n[3])+eigb*(qn  *m[3]+H   *a2*n[3]); \
    RDL[4][4]   =  eiga*a1*(H   *m[4]+qn  *n[4])+eigb*(qn  *m[4]+H   *a2*n[4]); \
    \
    RDL[0][0]   += eig[0]; \
    RDL[1][1]   += eig[0]; \
    RDL[2][2]   += eig[0]; \
    RDL[3][3]   += eig[0]; \
    RDL[4][4]   += eig[0]; \
    \
    spra=  spra*e[3]; \
}

#define reigl_vis(u,mul,mut,e,tmpa,spra) \
{ \
	swFloat un,ke,mu,mupr,kq,rtmp; \
	int k; \
    mupr    =  mul/prl+mut/prt; \
    un      =  e[0]*u[1]+e[1]*u[2]+e[2]*u[3]; \
    ke      =  0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3]); \
    mu      =  mul+mut; \
    kq      =  cp*mupr; \
    rtmp    =  kq/(mu*rr); \
    \
    tmpa[0][0] =  0.0; \
    tmpa[1][0] = -u[1]-un*e[0]/3.0; \
    tmpa[2][0] = -u[2]-un*e[1]/3.0; \
    tmpa[3][0] = -u[3]-un*e[2]/3.0; \
    tmpa[4][0] = -rtmp*u[4]/u[0]-un*un/3.0-ke*2.0+gk1*rtmp*ke; \
    tmpa[0][1] =  0.0; \
    tmpa[1][1] =  e[0]*e[0]/3.0+1.0; \
    tmpa[2][1] =  e[0]*e[1]/3.0; \
    tmpa[3][1] =  e[0]*e[2]/3.0; \
    tmpa[4][1] =  u[1]+un*e[0]/3.0-gk1*rtmp*u[1]; \
    tmpa[0][2] =  0.0; \
    tmpa[1][2] =  e[0]*e[1]/3.0; \
    tmpa[2][2] =  e[1]*e[1]/3.0+1.0; \
    tmpa[3][2] =  e[1]*e[2]/3.0; \
    tmpa[4][2] =  u[2]+un*e[1]/3.0-gk1*rtmp*u[2]; \
    tmpa[0][3] =  0.0; \
    tmpa[1][3] =  e[0]*e[2]/3.0; \
    tmpa[2][3] =  e[1]*e[2]/3.0; \
    tmpa[3][3] =  e[2]*e[2]/3.0+1.0; \
    tmpa[4][3] =  u[3]+un*e[2]/3.0-gk1*rtmp*u[3]; \
    tmpa[0][4] =  0.0; \
    tmpa[1][4] =  0.0; \
    tmpa[2][4] =  0.0; \
    tmpa[3][4] =  0.0; \
    tmpa[4][4] =  gk1*rtmp; \
    rtmp       =  e[3]*mu*e[4]/u[0]; \
    for(j=0;j<5;j++) \
    { \
    	for(k=0;k<5;k++) \
    	{ \
    		tmpa[j][k] = tmpa[j][k]*rtmp; \
    	} \
    } \
    spra=  e[3]*e[4]*gk*mupr/u[0]; \
    spra=  tmpa[4][0]; \
}

//define 2 function pointers
define_e2v_FunPtr(imp_lhs_imp)
{
	//selfConn computation
	int iedge,j,idxL,idxR;
	swFloat* fv_u      = accessArray(data->vertexData, 0);
	swFloat* sec_bct   = accessArray(data->vertexData, 1);
	swFloat* fv_mu     = accessArray(data->vertexData, 2);
	swInt dims_mu      = getArrayDims(data->vertexData, 2);
	swFloat* fv_uc     = accessArray(data->vertexData, 3);

	swFloat* n_vg      = accessArray(data->frontEdgeData, 0);
	swFloat* transform = accessArray(data->frontEdgeData, 1);

	swFloat* para      = accessArray(paraData, 0);

	swInt edgeNum = getArraySize(data->frontEdgeData);

	swFloat *nn, *vg, *uL, *uR;
	swFloat JL[5][5], JR[5][5], spra, eig[3], um[5], mul, mut, JV[5][5], sprv;
	swFloat *JV_ptr = &JV[0][0];
	swFloat *JL_ptr = &JL[0][0];
	swFloat *JR_ptr = &JR[0][0];
	swFloat gk  = para[2];
	swFloat gk1 = para[3];
	swFloat LHS_efix = para[4];
	swFloat prl = para[5];
	swFloat prt = para[6];
	swFloat rr  = para[7];
	swFloat cp  = para[8];
	for(iedge=0;iedge<edgeNum;iedge++)
	{
		if(transform[iedge]==0)
		{
			idxL = startVertices[iedge];
			idxR = endVertices[iedge];
		} else
		{
			idxL = endVertices[iedge];
			idxR = startVertices[iedge];
		}

		nn = &n_vg[iedge*8+0];
		vg = &n_vg[iedge*8+5];
		uL = &fv_u[idxL*5+0];
		uR = &fv_u[idxR*5+0];

		reigl_new(uL, vg, nn,  1.0, 1.0, LHS_efix, JL, spra, eig)
		reigl_new(uR, vg, nn, -1.0, 1.0, LHS_efix, JR, spra, eig)

        if(para[0]==1)
        {
            um[0]  =  0.5*(uL[0]+uR[0]);
            um[1]  =  0.5*(uL[1]+uR[1]);
            um[2]  =  0.5*(uL[2]+uR[2]);
            um[3]  =  0.5*(uL[3]+uR[3]);
            um[4]  =  0.5*(uL[4]+uR[4]);
            mul    =  0.5*(fv_mu[idxL*dims_mu+0]+fv_mu[idxR*dims_mu+0]);
            if(para[1]==1)  mut =  0.5*(fv_mu[idxL*dims_mu+1]+fv_mu[idxR*dims_mu+1]);
            reigl_vis(um, mul, mut, nn, JV, sprv);
            for(j=0;j<25;j++)
            {
            	JL_ptr[j] += JV_ptr[j];
            	JR_ptr[j] -= JV_ptr[j];
            }
        }

		fv_uc[idxL*5+0] += JL[0][0]+JR[0][1];
		fv_uc[idxL*5+1] += JL[1][0]+JR[1][1];
		fv_uc[idxL*5+2] += JL[2][0]+JR[2][1];
		fv_uc[idxL*5+3] += JL[3][0]+JR[3][1];
		fv_uc[idxL*5+4] += JL[4][0]+JR[4][1];
	}
}

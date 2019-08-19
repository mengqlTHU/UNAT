#include "swMacro.h"
#include "iterator.h"
#include "math.h"

#define rhs_conv_roe(calD, addD) \
{ \
	/* 273 flops */ \
	swFloat duc[5],u[5],eig[3],e[3],ad[5],tv[5]; \
	swFloat aL,aR,hL,hR,unL,unR,d,h,a,ke,un,vn,area,uan,eiga,eigb,c1,c2; \
	swFloat tmp1,tmp2,tmp3; \
	/* 31 flops */ \
	tmp1 = 0.5*(uLf[1]*uLf[1]+uLf[2]*uLf[2]+uLf[3]*uLf[3]); \
	tmp2 = 0.5*(uRf[1]*uRf[1]+uRf[2]*uRf[2]+uRf[3]*uRf[3]); \
	duc[0] = uRf[0]-uLf[0]; \
	duc[1] = uRf[0]*uRf[1]-uLf[0]*uLf[1]; \
	duc[2] = uRf[0]*uRf[2]-uLf[0]*uLf[2]; \
	duc[3] = uRf[0]*uRf[3]-uLf[0]*uLf[3]; \
	duc[4] = (uRf[4]-uLf[4])/0.4+uRf[0]*tmp2-uLf[0]*tmp1; \
	/* 34 flops*/ \
	aL = sqrt(1.4*uLf[4]/uLf[0]); \
	hL = aL*aL/0.4+tmp1; \
	aR = sqrt(1.4*uRf[4]/uRf[0]); \
	hR = aR*aR/0.4+tmp2; \
	/* 58 flops*/ \
	u[0] = sqrt(uLf[0]*uRf[0]); \
	d    = sqrt(uRf[0]/uLf[0]); \
	tmp1 = 1.0/(1.0+d); \
	u[1] = (uLf[1]+d*uRf[1])*tmp1; \
	u[2] = (uLf[2]+d*uRf[2])*tmp1; \
	u[3] = (uLf[3]+d*uRf[3])*tmp1; \
	h    = (hL+d*hR)*tmp1; \
	ke   = 0.5*(u[1]*u[1]+u[2]*u[2]+u[3]*u[3]); \
	u[4] = (h-ke)*0.4*u[0]/1.4; \
	a    = sqrt(fabs(1.4*u[4]/u[0])); \
	/* 17 flops*/\
	e[0] = fac_1d[0]; \
	e[1] = fac_1d[1]; \
	e[2] = fac_1d[2]; \
	area = fac_1d[3]; \
	vn   = vgf[0]*e[0]+vgf[1]*e[1]+vgf[2]*e[2]; \
	unL = e[0]*uLf[1]+e[1]*uLf[2]+e[2]*uLf[3]-vn; \
	unR = e[0]*uRf[1]+e[1]*uRf[2]+e[2]*uRf[3]-vn; \
	/* 27 flops */\
	if(1) \
	{ \
		tmp1  = uLf[0]*unL; \
		tmp2  = uRf[0]*unR; \
		tv[0] = tmp1+tmp2; \
		tv[1] = tmp1*uLf[1]+tmp2*uRf[1]+e[0]*(uLf[4]+uRf[4]); \
		tv[2] = tmp1*uLf[2]+tmp2*uRf[2]+e[1]*(uLf[4]+uRf[4]); \
		tv[3] = tmp1*uLf[3]+tmp2*uRf[3]+e[2]*(uLf[4]+uRf[4]); \
		tv[4] = tmp1*hL    +tmp2*hR    +vn  *(uLf[4]+uRf[4]); \
	} else \
	{ \
	} \
	/* 35 flops */\
	uan    = u[1]*e[0]+u[2]*e[1]+u[3]*e[2]; \
	un     = uan-vn; \
	eig[0] = fabs(un); \
	eig[1] = fabs(un+a); \
	eig[2] = fabs(un-a); \
	tmp1   = 2.0*fabs(unR+aR-unL-aL); \
	if(eig[1]<tmp1) eig[1] = 0.5*(eig[1]*eig[1]/tmp1+tmp1); \
	tmp1   = 2.0*fabs(unR-aR-unL+aL); \
	if(eig[2]<tmp1) eig[2] = 0.5*(eig[2]*eig[2]/tmp1+tmp1); \
	/* 16 flops */\
	eiga = 0.5*(eig[1]+eig[2])-eig[0]; \
	eigb = 0.5*(eig[1]-eig[2]); \
	tmp1 = -uan*duc[0]+e[0]*duc[1]+e[1]*duc[2]+e[2]*duc[3]; \
	tmp2 = ke*duc[0]-u[1]*duc[1]-u[2]*duc[2]-u[3]*duc[3]+duc[4]; \
	tmp2 = tmp2*0.4; \
	/* 55 flops */\
	a  = 1.0/a; \
	c1 = (eiga*tmp2*a+eigb*tmp1)*a; \
	c2 = eigb*tmp2*a+eiga*tmp1; \
	ad[0] = eig[0]*duc[0]; \
	ad[1] = eig[0]*duc[1]; \
	ad[2] = eig[0]*duc[2]; \
	ad[3] = eig[0]*duc[3]; \
	ad[4] = eig[0]*duc[4]; \
	\
	ad[0] = ad[0]+c1; \
	ad[1] = ad[1]+c1*u[1]+c2*e[0]; \
	ad[2] = ad[2]+c1*u[2]+c2*e[1]; \
	ad[3] = ad[3]+c1*u[3]+c2*e[2]; \
	ad[4] = ad[4]+c1*h   +c2*uan; \
	if((addD)==1) \
	{ \
		rhsl[0] = 0.5*area*(tv[0]-ad[0]); \
		rhsl[1] = 0.5*area*(tv[1]-ad[1]); \
		rhsl[2] = 0.5*area*(tv[2]-ad[2]); \
		rhsl[3] = 0.5*area*(tv[3]-ad[3]); \
		rhsl[4] = 0.5*area*(tv[4]-ad[4]); \
	} else \
	{ \
		rhsl[0] = 0.5*area*tv[0]; \
		rhsl[1] = 0.5*area*tv[1]; \
		rhsl[2] = 0.5*area*tv[2]; \
		rhsl[3] = 0.5*area*tv[3]; \
		rhsl[4] = 0.5*area*tv[4]; \
		\
		rhsD[0] = 0.5*area*ad[0]; \
		rhsD[1] = 0.5*area*ad[1]; \
		rhsD[2] = 0.5*area*ad[2]; \
		rhsD[3] = 0.5*area*ad[3]; \
		rhsD[4] = 0.5*area*ad[4]; \
	} \
}


//define 2 function pointers
define_e2v_FunPtr(upw_conv)
{
	//selfConn computation
	int i,j;
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
	swFloat* order     = accessArray(data->vertexData, 5);
	swFloat* is_int    = accessArray(data->vertexData, 6);
	
	//frontEdge computation
	swFloat* n_vg	      = accessArray(data->frontEdgeData, 0);
	swInt dims_vg         = getArrayDims(data->frontEdgeData, 0);
	swFloat* transform    = accessArray(data->frontEdgeData, 1);
	swInt dims_trans      = getArrayDims(data->frontEdgeData, 1);
	swFloat* mortar_cen   = accessArray(data->frontEdgeData, 2);
	swInt dims_mortar_cen = getArrayDims(data->frontEdgeData, 2);

	swFloat* para = accessArray(paraData, 0);


	swInt edgeNumber = getArraySize(data->frontEdgeData);
	swInt iedge,idxL,idxR;
	swFloat fac_1d[4], vgf[3], uLf[5],uRf[5], dL[3], dR[3], duL[5], duR[5];
	swFloat limiter_L=1.0, limiter_R=1.0;
	swFloat rhsl[5],rhsD[5];
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		if(transform[iedge] == 0)	
		{
			idxL = startVertices[iedge];
			idxR = endVertices[iedge];
		} else if(transform[iedge] == 1)
		{
			idxL = endVertices[iedge];
			idxR = startVertices[iedge];
			// printf("%d\n",iedge);
		} else
		{
			// idxL = endVertices[iedge];
			// idxR = startVertices[iedge];
			// printf("%d\n",(int)transform[iedge]);
			printf("error\n");
		}
		fac_1d[0] = n_vg[iedge*dims_vg+0];
		fac_1d[1] = n_vg[iedge*dims_vg+1];
		fac_1d[2] = n_vg[iedge*dims_vg+2];
		fac_1d[3] = n_vg[iedge*dims_vg+3];
		vgf[0]    = n_vg[iedge*dims_vg+5];
		vgf[1]    = n_vg[iedge*dims_vg+6];
		vgf[2]    = n_vg[iedge*dims_vg+7];

		uLf[0] = fv_u[idxL*dims_u+0];
		uLf[1] = fv_u[idxL*dims_u+1];
		uLf[2] = fv_u[idxL*dims_u+2];
		uLf[3] = fv_u[idxL*dims_u+3];
		uLf[4] = fv_u[idxL*dims_u+4];

		uRf[0] = fv_u[idxR*dims_u+0];
		uRf[1] = fv_u[idxR*dims_u+1];
		uRf[2] = fv_u[idxR*dims_u+2];
		uRf[3] = fv_u[idxR*dims_u+3];
		uRf[4] = fv_u[idxR*dims_u+4];

		if(order[idxL] > 1)
		{

			dL[0] = mortar_cen[iedge*dims_mortar_cen+0]-sec_cen[idxL*dims_sec_cen+0];
			dL[1] = mortar_cen[iedge*dims_mortar_cen+1]-sec_cen[idxL*dims_sec_cen+1];
			dL[2] = mortar_cen[iedge*dims_mortar_cen+2]-sec_cen[idxL*dims_sec_cen+2];

			dR[0] = mortar_cen[iedge*dims_mortar_cen+0]-sec_cen[idxR*dims_sec_cen+0];
			dR[1] = mortar_cen[iedge*dims_mortar_cen+1]-sec_cen[idxR*dims_sec_cen+1];
			dR[2] = mortar_cen[iedge*dims_mortar_cen+2]-sec_cen[idxR*dims_sec_cen+2];

			// 190 flops
			if(para[1]==1)
			// if(0)
			{
				limiter_L = fv_gra[idxL*dims_gra+18];
				limiter_R = fv_gra[idxR*dims_gra+18];
			}
			for(j=0;j<5;j++)
			{
				duL[j] = fv_gra[idxL*dims_gra+3*j+0]*dL[0]
					+ fv_gra[idxL*dims_gra+3*j+1]*dL[1]
					+ fv_gra[idxL*dims_gra+3*j+2]*dL[2];
				duR[j] = fv_gra[idxR*dims_gra+3*j+0]*dR[0]
					+ fv_gra[idxR*dims_gra+3*j+1]*dR[1]
					+ fv_gra[idxR*dims_gra+3*j+2]*dR[2];
				if(order[idxL] == 2)
				// if(1)
				{
					uLf[j] += limiter_L*duL[j];
					uRf[j] += limiter_R*duR[j];
				} else
				{

				}
			}
			if((uLf[0] <= 0 || uLf[4] <= 0) || (uRf[0] <= 0 || uRf[4] <= 0))
			{
				for (j = 0; j < 5; j++)
				{
					uLf[j] = fv_u[idxL*dims_u+j];
					uRf[j] = fv_u[idxR*dims_u+j];
				}
			}
		}
		
		// 273 flops
		if(para[0]==2)
		// if(1)
		{
			rhs_conv_roe(para[2], para[3]);
			// rhs_conv_roe(1, 0);
		}

		// 15 flops
		for(j=0;j<5;j++)
		{
			fv_rhs[idxL*dims_rhs+j] += rhsl[j];
		}
		if(is_int[idxR]==1)
		{
			for(j=0;j<5;j++)
			{
				fv_rhs[idxR*dims_rhs+j] -= rhsl[j];
			}
		}
		// if(0)
		if(para[2]==1 && para[3]==-1)
		{
			for(j=0;j<5;j++)
			{
				fv_duc[idxL*dims_duc+j] += rhsD[j];
			}
			if(is_int[idxR]==1)
			{
				for(j=0;j<5;j++)
				{
					fv_duc[idxR*dims_duc+j] -= rhsD[j];
				}
			}
		}
// if(iedge==0) printf("%f,%d,%d,%f,%f\n",transform[iedge],idxL,idxR,uLf[1],uRf[1]);
	}
}



/*
* @Author: Gao Fei
* @Date:   2018-12-07 15:22:32
* @Last Modified by:   Gao Fei
* @Last Modified time: 2018-12-07 16:44:10
*/

#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <vector>
#include <cmath>
#include <assert.h>
#include "swMacro.h"
#include "topology.H"
#include "iterator.H"
#include "iterator.h"
#include "directSegmentIterator.H"
#include "fluxMacros.h"
#include "funMacros.h"

#define NONZERONUM 1516800
#define MTXORDRNUM 512000

int* readFile(char* name);
void checkResults(swFloat* array1, swFloat* array2, swInt count);

void main()
{
	// read basic
	char owner[] = "owner_512000";
	char neighbor[] = "neighbour_512000";
	swInt* own = readFile(owner);
	swInt* nei = readFile(neighbor);

	swInt faceNum = NONZERONUM;
	swInt cellNum = MTXORDRNUM;

	swInt* cellWeights = new swInt[cellNum];
	swInt* edgeWeights = new swInt[faceNum];
	for (swInt i = 0; i < cellNum; ++i)
		cellWeights[i] = 0;
    for (swInt i = 0; i < faceNum; ++i)
    {
        edgeWeights[i] = 1;
        cellWeights[ own[i] ] ++;
        cellWeights[ nei[i] ] ++;
    }

	// construct topology
	Topology* topo = Topology::constructFromEdge(own,nei,faceNum);
	if( topo->getVertexNumber() != cellNum) printError("wrong vertex number!");

	// define iterator
	DirectSegmentIterator iterator(*topo, &cellWeights[0], &edgeWeights[0]);

	checkDSI(iterator);

	snGrad(t);
}

void snGrad(Topology* topo, DirectSegmentIterator& iterator, swInt nCpnt)//swInt *own, swInt *nei,
{
	// basic
	swInt faceNum = topo->getEdgeNumber();
	swInt cellNum = topo->getVertexNumber();

	// declare data arrays
	swFloat* ssf = new swFloat[ faceNum*nCpnt ];
	swFloat* ssfRight
				 = new swFloat[ faceNum*nCpnt ];
	swFloat* deltaCoeffs
				 = new swFloat[ faceNum ];
	swFloat* vf  = new swFloat[ cellNum*nCpnt ];

	// initiliza data arrays
	for (swInt i = 0; i < cellNum; ++i)
		for (swInt j = 0; j < nCpnt; ++j) vf[ i*nCpnt + j ] = 1.0;
	for (swInt i = 0; i < faceNum; ++i) deltaCoeffs[i] = 1.0;

	// MPE computation for verification and comparation
	for (swInt i = 0; i < faceNum; ++i)
	{
		for (swInt j = 0; j < nCpnt; ++j)
			ssfRight[ i*nCpnt + j ]
				= deltaCoeffs[i]
				 * ( vf[ nei[i]*nCpnt + j ] - vf[ own[i]*nCpnt + j ] );
	}


	// CPE computation
	// single constructor for arrays
	Arrays backEdgeData, frontEdgeData, selfConnData, vertexData;

	// needed but not used
	constructEmptyArray( backEdgeData );
	constructEmptyArray( selfConnData );

	// face datas
	constructSingleArray(frontEdgeData, nCpnt, faceNum, COPYOUT, ssf);?????????????????????????????????????????
	addSingleArray(frontEdgeData, 1, faceNum, COPYIN, deltaCoeffs);

	// cell datas
	constructSingleArray(vertexData, nCpnt, cellNum, COPYIN, vf);

    // compute results with iterator
	iterator.edge2VertexIteration( &backEdgeData, &frontEdgeData,
		&selfConnData, &vertexData, fvcSnGrad_host, fvcSnGrad_slave);


	// checkResults
	checkResults(ssf, ssfRight, nCpnt*faceNum);

	// delete
	delete [] ssf;
	delete [] ssfRight;
	delete [] deltaCoeffs;
	delete [] vf;
}

void surfaceIntegrate(Topology& t, swInt *own, swInt *nei, swInt nCpnt)
{
	// basic
	swInt faceNum = t.getEdgeNumber();
	swInt cellNum = t.getVertexNumber();

	// declare data arrays
	swFloat* ivf = new swFloat[ nCpnt*cellNum ];
	swFloat* ivfRight= new swFloat[ nCpnt*cellNum ];
	swFloat* issf = new swFloat[ nCpnt*faceNum ];

	// initilize data arrays
	for (swInt i = 0; i < cellNum; ++i)
	{
		for (swInt j = 0; j < nCpnt; ++j)
		{
			ivf[ i*nCpnt + j ] = 0.0;
			ivfRight[ i*nCpnt + j ] = 0.0;
		}
	}
	for (swInt i = 0; i < faceNum; ++i)
	{
		for (swInt j = 0; j < nCpnt; ++j) issf[ i*nCpnt + j ] = (faceNum%7)*1.0;
	}

	// MPE computation for verification and comparation
	for (swInt i = 0; i < faceNum; ++i)
	{
		for (swInt j = 0; j < nCpnt; ++j)
		{
			ivfRight[ own[i]*nCpnt + j ] += issf[ i*nCpnt + j ];
			ivfRight[ nei[i]*nCpnt + j ] -= issf[ i*nCpnt + j ];
		}
	}

	// CPE computation
	// single constructor for arrays
	Arrays backEdgeData, frontEdgeData, selfConnData, vertexData;

	// needed but not used
	constructEmptyArray( backEdgeData );
	constructEmptyArray( selfConnData );

	// face datas
	constructSingleArray(frontEdgeData, nCpnt, faceNum, COPYOUT, ssf);?????????????????????????????????????????
	addSingleArray(frontEdgeData, 1, faceNum, COPYIN, deltaCoeffs);

	// cell datas
	constructSingleArray(vertexData, nCpnt, cellNum, COPYIN, vf);

    // compute results with iterator
	iterator.edge2VertexIteration( &backEdgeData, &frontEdgeData,
		&selfConnData, &vertexData, fvcSnGrad_host, fvcSnGrad_slave);

	// checkResults
	checkResults(ivf, ivfRight, nCpnt*cellNum);

	// delete

}

void negSumDiag(Topology& t, swInt *own, swInt *nei)
{
	// basic
	swInt faceNum = t.getEdgeNumber();
	swInt cellNum = t.getVertexNumber();

	// declare data arrays
	swFloat* Diag = new swFloat[ cellNum ];
	swFloat* DiagRight = new swFloat[ cellNum ];
	swFloat* Lower = new swFloat[ faceNum ];
	swFloat* Upper = new swFloat[ faceNum ];

	// initilize data arrays
	for (swInt i = 0; i < faceNum; ++i)
	{
		Lower[i] = i%7*1.0;
		Upper[i] = i%13*1.0;
	}
	for (swInt i = 0; i < cellNum; ++i)
	{
		Diag[i] = 0.0;
	}

	// MPE computation for verification and comparation
	for (swInt i = 0; i < faceNum; ++i)
	{
		DiagRight[ own[i] ] -= Lower[i];
		DiagRight[ nei[i] ] -= Upper[i];
	}

	// CPE computation

	// checkResults
	checkResults(Diag, DiagRight, cellNum);

	// delete

}

void calcLimiter(Topology& t, swInt *own, swInt *nei, swInt phiDim)
{
	// basic
	swInt faceNum = t.getEdgeNumber();
	swInt cellNum = t.getVertexNumber();

	// declare data arrays
	swFloat* limiter = new swFloat[ faceNum ];
	swFloat* limiterRight = new swFloat[ faceNum ];
	swFloat* CDweights = new swFloat[ faceNum ];
	swFloat* faceFlux_ = new swFloat[ faceNum ];
	if (phiDim == 1) {
		// scalar
		swFloat* lPhi = new swFloat[cellNum];
		// vector
		swFloat* gradc = new swFloat[spaceDim * cellNum];
	}
	else if (phiDim == 3) {
		//vector
		swFloat* lPhi = new swFloat[spaceDim * cellNum];
		//tensor
		swFloat* gradc = new swFloat[spaceDim * spaceDim * cellNum];
	}
	swFloat* C  = new swFloat[spaceDim * cellNum];

	// initilize data arrays
	for (swInt i = 0; i < faceNum; ++i)
	{
		CDweights[i] = ;
		faceFlux_[i] = ;
		C[i] = ;
	}
	for (swInt i = 0; i < cellNum; ++i)
	{
		if ( phiDim == 1 )
		{
			lPhi[i] = ;
			for (swInt j = 0; j<spaceDim)gradc[i*spaceDim + j] = ;
		}
	}

	// MPE computation for verification and comparation
	for (swInt i < faceNum; ++i)
	{
		limiterRight = limiter(
			CDweights[i],
			faceFlux_[i],
			lPhi[own[i]],
			lPhi[nei[i]],
			gradc + own[i]*spaceDim,
			gradc + nei[i]*spaceDim,
			C + own[i]*spaceDim,
			C + nei[i]*spaceDim,
			phiDim
			);
	}

	// CPE computation

	// checkResults


	// delete
}

swFloat limiter(
	swFloat CDweights,
	swFloat faceFlux,
	swFloat phiP,
	swFloat phiN,
	swFloat* gradPPtr,
	swFloat* gradNPtr,
	swFloat* ownC,
	swFloat* neiC,
	swInt phiDim
	)
{
	// calculate slope rate
	swFloat d[spaceDim] = {neiC[0]-ownC[0], neiC[1]-ownC[1], neiC[2]-ownC[2]};
	swFloat gradf = phiN - phiP;
	swFloat gradcf;
	if (faceFlux > 0)
    {
        gradcf = d[0]*gradcP[0] + d[1]*gradcP[1] + d[2]*gradcP[2];
    }
    else
    {
        gradcf = d[0]*gradcN[0] + d[1]*gradcN[1] + d[2]*gradcN[2];
    }

    if (fabs(gradcf) >= 1000*fabs(gradf))
    {
        r = 2*1000*sign(gradcf)*sign(gradf) - 1;
    }
    else
    {
        r = 2*(gradcf/gradf) - 1;
    }

    swFloat k = 1.0;
    swFloat twoByk_ = 2.0/max(k, SMALL);

    return max(min(twoByk_*r, 1.0), 0.0);

	// for limitedLinear scheme
}

swFloat max(swFloat a, swFloat b)
{
	if (a>b) return a;
	else return b;
}

swFloat min(swFloat a, swFloat b)
{
	if (a<b) return a;
	else return b;
}

swFloat sign(swFloat a)
{
	if(a>0.0) return 1.0;
	else return -1.0;
}

void linearUpwind(Topology& t, swInt *own, swInt *nei, swInt phiDim)
{
	// basic
	swInt faceNum = t.getEdgeNumber();
	swInt cellNum = t.getVertexNumber();

	// declare data arrays
	if ( phiDim == 3 )
	{
		swFloat* sfCorr = new swFloat[spaceDim*faceNum];
		swFloat* sfCorrRight = swFloat[spaceDim*faceNum];
		swFloat* gradVf = new swFloat[spaceDim * spaceDim * cellNum];
	}
	else {
		// print error
	}

	swFloat* C  = new swFloat[spaceDim * cellNum];
	swFloat* Cf = new swFloat[spaceDim * faceNum];
	swFloat* faceFlux = new swFloat[faceNum];

	// initilize
	for (swInt i = 0; i < faceNum; ++i)
	{
		faceFlux[i] = ;
		for (swInt j = 0; j < spaceDim; ++j)
		{
			Cf[i*spaceDim + j] = ;
			sfCorrRight[ i*spaceDim+j ] = 0.0;
			sfCorrt[ i*spaceDim+j ] = 0.0;
		}
	}

	for (swInt i = 0; i < cellNum; ++i)
	{
		for (swInt j = 0; j < spaceDim; ++j)
		{
			C[i*spaceDim + j] = ;
			for (swInt k = 0; k < spaceDim; ++k)
			{
				gradVf[ (i*spaceDim + j)*spaceDim + k] = ;
			}
		}
	}

	// MPE computation for verification and comparation
	for (swInt i = 0; i < faceNum; ++i)
	{
		swInt celli = (faceFlux[i]) ? own[i] : nei[i];
		for (swInt j = 0; j < spaceDim; ++j)
		{
			for (swInt k = 0; k < spaceDim; ++k)
				sfCorrRight[ i*spaceDim+j ]
					= ( Cf[ i*spaceDim+k ] - C[ celli*spaceDim+k ] )
					 *  gradVf[ (i*spaceDim + j)*spaceDim + k];
		}
	}

	// CPE computation

	// checkResults

	// delete

}

void cellLimited_calcMaxMin(Topology& t, swInt *own, swInt *nei, swInt nCpnt)
{
	// basic
	swInt faceNum = t.getEdgeNumber();
	swInt cellNum = t.getVertexNumber();

	// declare data arrays
	swFloat* maxVsf = new swFloat[ cellNum*nCpnt ];
	swFloat* maxVsfRight = new swFloat[ cellNum*nCpnt ];
	swFloat* minVsf = new swFloat[ cellNum*nCpnt ];
	swFloat* minVsfRight = new swFloat[ cellNum*nCpnt ];
	swFloat* vsf = new swFloat[ cellNum*nCpnt ];

	// initilize data arrays
	for (swInt i = 0; i < cellNum; ++i)
	{
		for (swInt j = 0; j < nCpnt; ++j)
		{
			vsf[i*spaceDim + j] = ((i*j)%7)*1.0;
			maxVsf[i*spaceDim + j] = -8.0;
			maxVsfRight[i*spaceDim + j] = -8.0;
			minVsf[i*spaceDim + j] = 8.0;
			minVsfRight[i*spaceDim + j] = 8.0;
		}
	}

	// MPE computation for verification and comparation
	for (swInt i = 0; i < faceNum; ++i)
	{
		swFloat* vsfOwnPtr = vsf + own[i]*nCpnt;
		swFloat* vsfNeiPtr = vsf + nei[i]*nCpnt;
		swFloat* maxVsfOwnPtr = maxVsfRight + own[i]*nCpnt;
		swFloat* maxVsfNeiPtr = maxVsfRight + nei[i]*nCpnt;
		swFloat* minVsfOwnPtr = minVsfRight + own[i]*nCpnt;
		swFloat* minVsfNeiPtr = minVsfRight + nei[i]*nCpnt;

		max_v(maxVsfOwnPtr, vsfNeiPtr, nCpnt);
		min_v(minVsfOwnPtr, vsfNeiPtr, nCpnt);
		max_v(maxVsfNeiPtr, vsfOwnPtr, nCpnt);
		min_v(minVsfOwnPtr, vsfOwnPtr, nCpnt);
	}

	// CPE computation

	// checkResults
	checkResults(maxVsf, maxVsfRight, cellNum*nCpnt);
	checkResults(minVsf, minVsfRight, cellNum*nCpnt);

	// delete
}

void max_v(swFloat* real1, swFloat* real2, swInt nCpnt)
{
	for (swInt i = 0; i < nCpnt; ++i)
	{
		real1[i]= (real1[i]>real2[i]) ? real1[i] :  real2[i];
	}
}
void min_v(swFloat* real1, swFloat* real2, swInt nCpnt)
{
	for (swInt i = 0; i < nCpnt; ++i)
	{
		real1[i]= (real1[i]<real2[i]) ? real1[i] :  real2[i];
	}
}

void cellLimited_limitFace(Topology& t, swInt *own, swInt *nei, swInt nCpnt)
{
	// basic
	swInt faceNum = t.getEdgeNumber();
	swInt cellNum = t.getVertexNumber();

	// declare data arrays
	swFloat* limiter = new swFloat[ cellNum*nCpnt ];
	swFloat* limiterRight = new swFloat[ cellNum*nCpnt ];
	swFloat* C  = new swFloat[spaceDim * cellNum];
	swFloat* Cf = new swFloat[spaceDim * faceNum];
	swFloat* maxVsf = new swFloat[ cellNum*nCpnt ];
	swFloat* minVsf = new swFloat[ cellNum*nCpnt ];
	swFloat* g = new swFloat[ spaceDim*cellNum*nCpnt ];

	// initiliza data arrays
	for (swInt i = 0; i < cellNum; ++i)
		for (swInt j = 0; j < nCpnt; ++j)
		{
			limiter[ i*nCpnt + j ] = 2.0;
			limiterRight[ i*nCpnt + j ] = 2.0;
			maxVsf[ i*nCpnt + j ] =  100.0;
			minVsf[ i*nCpnt + j ] = -100.0;
		}
	for (swInt i = 0; i < faceNum; ++i)
		for (swInt j = 0; j < spaceDim; ++j)
			Cf[i*spaceDim + j] = ;
	for (swInt i = 0; i < cellNum; ++i)
		for (swInt j = 0; j < spaceDim; ++j)
			C[i*spaceDim + j] = ;

	for (swInt i = 0; i < cellNum; ++i)
		for (swInt j = 0; j < nCpnt; ++j)
			for (swInt k = 0; k < spaceDim; ++k)
				g[(i*spaceDim + j)*spaceDim + k]=1.0+k*1.0;

	// MPE computation for verification and comparation
	for (swInt i = 0; i < faceNum; ++i)
	{
		for (swInt j = 0; j < nCpnt; ++j)
		{
			ssfRight[ i*nCpnt + j ];
			for (swInt k = 0; k < spaceDim; ++k)
				g[(i*spaceDim + j)*spaceDim + k]
		}
	}

	// CPE computation

	// checkResults
	checkResults(ssf, ssfRight, nCpnt*faceNum);

	// delete
}

void checkResults(swFloat* array1, swFloat* array2, swInt count)
{
	cout<<"check result..."<<endl;
	for(int i=0;i<count;i++)
	{
		if(fabs(array1[i]-array2[i])>EPS)
		{
			if(fabs(array1[i]-array2[i])/array1[i]>EPS)
			{
				cout<<"Error on index["<<i<<"], "
					<<array1[i]<<", "<<array2[i]<<endl;
				exit(1);
			}
		}
	}
	cout<<"The result is correct!"<<endl;
}

int* readFile(char* name){
	FILE *fp = fopen(name,"r");
	if(fp==NULL){
		return 0;
	}
	char c1[100],c2[100];
	int cellNum;
	while(1)
	{
		fgets(c1,100,fp);
		if(c1[0]=='('){
			cellNum = atoi(c2);
			break;
		}
		fgets(c2,100,fp);
		if(c2[0]=='('){
			cellNum = atoi(c1);
			break;
		}else{
		}
	}
	int *buff = (int*)malloc(sizeof(int)*NONZERONUM);
	int i=0;
	while(i<NONZERONUM)
	{
		fgets(c1,100,fp);
		buff[i] = atoi(c1);
		i++;
	}
	printf("InternalFaceNum = %d\n",i);
	return buff;
}

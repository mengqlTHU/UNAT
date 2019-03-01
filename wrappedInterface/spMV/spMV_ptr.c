#include "spMV.h"
#include "swMacro.h"
#include "iterator.h"

//define 2 function pointers
define_e2v_FunPtr(spMV)
{
	//selfConn computation
	swFloat* diag	= accessArray(selfConnData, 0);
	swFloat* b		= accessArray(vertexData, 0);
	swFloat* x		= accessArray(vertexData, 1);
	swInt iDim,dims;
	
	swInt vertexNum = getArraySize(selfConnData);
	dims = getArrayDims(selfConnData, 0);
	swInt ivertex;
	for( ivertex = 0; ivertex < vertexNum; ivertex++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			b[ivertex*dims+iDim]
				+= diag[ivertex*dims+iDim]*x[ivertex*dims+iDim];
		}
	}
	
	//frontEdge computation
	swFloat* upper	= accessArray(frontEdgeData, 0);
	swInt edgeNumber = getArraySize(frontEdgeData);
	dims = getArrayDims(frontEdgeData, 0);
	swInt iedge;
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
//if(startVertices[iedge]==105) printf("owner:%d,%d,%f,%f,%f\n",iedge,endVertices[iedge],b[startVertices[iedge]],upper[iedge],x[endVertices[iedge]]);
			b[startVertices[iedge]*dims+iDim]
				+= upper[iedge*dims+iDim]*x[endVertices[iedge]*dims+iDim];
		}
	}

	//backEdge computation
	swFloat* lower	= accessArray(backEdgeData, 0);
	dims = getArrayDims(backEdgeData, 0);
	edgeNumber = getArraySize( backEdgeData );
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
//if(endVertices[iedge]==105) printf("neighbor:%d,%d,%f,%f,%f\n",iedge,startVertices[iedge],b[endVertices[iedge]],lower[iedge],x[startVertices[iedge]]);
			b[endVertices[iedge]*dims+iDim]
				+= lower[iedge*dims+iDim]*x[startVertices[iedge]*dims+iDim];	
		}
	}
}



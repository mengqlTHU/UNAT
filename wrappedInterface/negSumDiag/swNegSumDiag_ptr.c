#include "swNegSumDiag.h"
#include <stdio.h>
#include "swMacro.h"
#include "iterator.h"

//define 2 function pointers
define_e2v_FunPtr(swNegSumDiag)
{
	swInt iDim,dims,iedge;
	swFloat *b;
	swInt edgeNumber = getArraySize(frontEdgeData);
	b = accessArray(vertexData,0);
	//frontEdge computation
	if(edgeNumber>0)
	{
	swFloat* upper	= accessArray(frontEdgeData, 0);
	dims = getArrayDims(frontEdgeData, 0);
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			b[startVertices[iedge]*dims+iDim]
				-= upper[iedge*dims+iDim];
		}
	}
	}

	//backEdge computation
	edgeNumber = getArraySize( backEdgeData );
	if(edgeNumber>0)
	{
	swFloat* lower	= accessArray(backEdgeData, 0);
	dims = getArrayDims(backEdgeData, 0);
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			b[endVertices[iedge]*dims+iDim]
				-= lower[iedge*dims+iDim];	
		}
	}
	}
}



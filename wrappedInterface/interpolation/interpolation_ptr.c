#include "interpolation.h"
#include <stdio.h>
#include "swMacro.h"
#include "iterator.h"

//define 2 function pointers
define_e2v_FunPtr(interpolation)
{
	swInt iDim,dims;
	swFloat* x = accessArray(vertexData,0);
	//frontEdge computation
	swFloat* ssf	= accessArray(frontEdgeData, 0);
	swFloat* coef	= accessArray(frontEdgeData, 1);
	swInt edgeNumber = getArraySize(frontEdgeData);
	dims = getArrayDims(frontEdgeData, 0);
	swInt iedge;
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
		for(iDim=0;iDim<dims;iDim++)
		{
			ssf[iedge*dims+iDim]
				= coef[iedge]
				* (x[endVertices[iedge]*dims+iDim]-x[startVertices[iedge]*dims+iDim]);
		}
	}
}



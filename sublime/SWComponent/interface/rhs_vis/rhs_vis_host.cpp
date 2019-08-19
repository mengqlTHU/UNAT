#include <stdio.h>
#include <stdlib.h>
#include "SW_iteratorInterface.h"
// #include "wrappedInterface.h"

#include "topology.hpp"
#include "iterator.hpp"
#include "multiLevelBlockIterator.hpp"
#include "rowSubsectionIterator.hpp"

extern "C" define_e2v_FunPtr(rhs_vis);
extern "C" define_e2v_slaveFunPtr(rhs_vis);
// extern void constructRSSIterator(swInt* rowAddr, swInt* colAddr,
	// swInt* cellWeights, swInt* faceWeights, swInt faceNum);
// extern void constructMLBIterator(swInt* rowAddr, swInt* colAddr,
	// swInt* cellWeights, swInt* faceWeights, swInt faceNum);
// 
extern "C"
{

 void rhs_vis_host_(swInt* rowAddr, swInt* colAddr,
 	swInt* faceNum, swInt* cellNum,
 	swFloat* transform, swFloat* n_vg, swFloat* u,
 	swFloat* sec_cen, swFloat* gra, swFloat* rhs, swFloat* duc,
 	swFloat* is_int, swFloat* mu, swFloat* turb, swFloat* t,
 	swFloat* addD, swFloat* is_tur_cal, swFloat* cp, swFloat* prL,
 	swFloat* prT, swFloat* is_KO, swFloat* is_RANS, swFloat* RANS_model)
 // void upw_conv_host_(swInt* faceNum)
 {
 	swInt face_num = *faceNum;
 	swInt cell_num = *cellNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 8, face_num, COPYIN, n_vg);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, transform);
 	constructSingleArray(vertexData, 5, cell_num, COPYOUT, rhs);
 	addSingleArray(vertexData, 5, cell_num, COPYIN, u);
 	addSingleArray(vertexData, 3, cell_num, COPYIN, sec_cen);
 	addSingleArray(vertexData, 19, cell_num, COPYIN, gra);
 	addSingleArray(vertexData, 5, cell_num, COPYOUT, duc);
 	addSingleArray(vertexData, 1, cell_num, COPYIN, is_int);
 	addSingleArray(vertexData, 1, cell_num, COPYIN, t);
 	int cellWeights = 39;
 	int faceWeights = 9;
 	if(*is_tur_cal==1) 
 	{
 		addSingleArray(vertexData, 2, cell_num, COPYIN, mu);
 		cellWeights += 2;
 	}
 	else
 	{
 		addSingleArray(vertexData, 1, cell_num, COPYIN, mu);
 		cellWeights += 1;
 	}
 	if(*is_RANS==1)
 	{
 		if(*RANS_model==1 || *RANS_model==2)
 		{
 			addSingleArray(vertexData, 1, cell_num, COPYIN, turb);
 			cellWeights += 1;
 		} else
 		{
 			addSingleArray(vertexData, 2, cell_num, COPYIN, turb);
 			cellWeights += 2;
 		}
 	}
 	constructEmptyArray(backEdgeData);
 	constructEmptyArray(selfConnData);

 	FieldData data;
 	data.backEdgeData = &backEdgeData;
 	data.frontEdgeData = &frontEdgeData;
 	data.vertexData = &vertexData;
 	data.selfConnData = &selfConnData;

 	coupledOperator opt;
 	opt.data      = &data;
 	opt.fun_host  = NULL;
 	opt.fun_slave = slave_rhs_vis;

 	Arrays paraData;
 	swFloat para[6] = {*addD, *is_tur_cal, *cp, *prL, *prT, *is_KO};
 	constructSingleArray(paraData, 1, 6, COPYIN, &para[0]);

 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo = UNAT::Topology::constructFromEdge(rowAddr, colAddr, face_num);
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, E2V, true);
	iter->edge2VertexIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
 }

}

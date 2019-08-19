#include <stdio.h>
#include <stdlib.h>
#include "swMacro.h"
#include "iterator.h"
#include "SW_iteratorInterface.h"

#include "topology.hpp"
#include "iterator.hpp"
#include "multiLevelBlockIterator.hpp"
#include "rowSubsectionIterator.hpp"
// #include "wrappedInterface.h"

extern "C" define_e2v_FunPtr(lhs_scalar);
extern "C" define_e2v_slaveFunPtr(lhs_scalar);

extern "C" define_e2v_FunPtr(lhs_scalar_zero);
extern "C" define_e2v_slaveFunPtr(lhs_scalar_zero);

extern "C" define_e2v_FunPtr(lhs_scalar_cal_dt);
extern "C" define_e2v_slaveFunPtr(lhs_scalar_cal_dt);

extern "C"
{

 void lhs_scalar_host_(swInt* rowAddr, swInt* colAddr,
 	swInt* faceNum, swInt* cellNum,
 	swFloat* transform, swFloat* n_vg, swFloat* u, swFloat* mu,
 	swFloat* LHS_s, swFloat* is_int, swFloat* is_vis_cal,
 	swFloat* is_tur_cal, swFloat* is_conv_only, swFloat* gk)
 // void upw_conv_host_(swInt* faceNum)
 {
 	// int i;
 	// printf("upw_conv_host: ");
 	// for(i=0;i<10;i++) printf("%f, ",transform[i]);
		// printf("\n");
 	// printf("upw_conv_host_\n");
 	swInt face_num = *faceNum;
 	swInt cell_num = *cellNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 8, face_num, COPYIN, n_vg);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, transform);
 	constructSingleArray(vertexData, 1, cell_num, COPYOUT, LHS_s);
 	addSingleArray(vertexData, 5, cell_num, COPYIN, u);
 	addSingleArray(vertexData, 1, cell_num, COPYIN, is_int);
  	int cellWeights = 7;
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
 	opt.fun_slave = slave_lhs_scalar;

 	Arrays paraData;
 	swFloat para[4] = {*is_vis_cal, *is_tur_cal, *is_conv_only, *gk};
 	constructSingleArray(paraData, 1, 4, COPYIN, &para[0]);

 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo = UNAT::Topology::constructFromEdge(rowAddr, colAddr, face_num);
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, E2V, true);
	iter->edge2VertexIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
 }
//void upw_conv_host_(){}

void lhs_scalar_zero_host_(swInt* faceNum, swFloat* LHS_s)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 1, face_num, COPYOUT, LHS_s);
 	constructEmptyArray(backEdgeData);
 	constructEmptyArray(selfConnData);
 	constructEmptyArray(vertexData);

 	FieldData data;
 	data.backEdgeData = &backEdgeData;
 	data.frontEdgeData = &frontEdgeData;
 	data.vertexData = &vertexData;
 	data.selfConnData = &selfConnData;

 	coupledOperator opt;
 	opt.data      = &data;
 	opt.fun_host  = NULL;
 	opt.fun_slave = slave_lhs_scalar_zero;

 	Arrays paraData;
 	constructEmptyArray(paraData);

 	int cellWeights = 1;
 	int faceWeights = 1;
 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo;
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, ARRAY, false);
	iter->arrayIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
}

void lhs_scalar_cal_dt_host_(swInt* faceNum, swFloat* LHS_s, swFloat* uc,
	swFloat* rhs, swFloat* is_int, swFloat* CFL, swFloat* eps)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 1, face_num, COPYINOUT, LHS_s);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, uc);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, rhs);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, is_int);
 	constructEmptyArray(backEdgeData);
 	constructEmptyArray(selfConnData);
 	constructEmptyArray(vertexData);

 	FieldData data;
 	data.backEdgeData = &backEdgeData;
 	data.frontEdgeData = &frontEdgeData;
 	data.vertexData = &vertexData;
 	data.selfConnData = &selfConnData;

 	coupledOperator opt;
 	opt.data      = &data;
 	opt.fun_host  = NULL;
 	opt.fun_slave = slave_lhs_scalar_cal_dt;

 	Arrays paraData;
 	swFloat para[2] = {*CFL, *eps};
 	constructSingleArray(paraData, 1, 2, COPYIN, &para[0]);

 	int cellWeights = 1;
 	int faceWeights = 12;
 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo;
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, ARRAY, false);
	iter->arrayIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
}

}



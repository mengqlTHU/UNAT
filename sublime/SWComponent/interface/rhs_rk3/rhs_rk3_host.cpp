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

/// Copy fv_uc to fv_uc0
extern "C" define_e2v_FunPtr(rhs_rk3_upt_rhs);
extern "C" define_e2v_slaveFunPtr(rhs_rk3_upt_rhs);

extern "C" define_e2v_FunPtr(rhs_rk3_zero);
extern "C" define_e2v_slaveFunPtr(rhs_rk3_zero);

extern "C"
{

void rhs_rk3_upt_rhs_host_(swInt* faceNum, swFloat* rv0, swFloat* duc,
	swFloat* rhs, swFloat* RK)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 5, face_num, COPYINOUT, rhs);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, duc);
 	addSingleArray(frontEdgeData, 5, face_num, COPYINOUT, rv0);
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
 	opt.fun_slave = slave_rhs_rk3_upt_rhs;

 	Arrays paraData;
 	swFloat para[1] = {*RK};
 	constructSingleArray(paraData, 1, 1, COPYIN, &para[0]);

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

void rhs_rk3_zero_host_(swInt* faceNum, swFloat* duc, swFloat* rhs)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 5, face_num, COPYOUT, duc);
 	addSingleArray(frontEdgeData, 5, face_num, COPYOUT, rhs);
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
 	opt.fun_slave = slave_rhs_rk3_zero;

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

}



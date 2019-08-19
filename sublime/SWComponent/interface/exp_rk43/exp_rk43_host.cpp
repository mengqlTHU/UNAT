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
extern "C" define_e2v_FunPtr(exp_rk43_copy);
extern "C" define_e2v_slaveFunPtr(exp_rk43_copy);

/// The 1st step of RK4 method
extern "C" define_e2v_FunPtr(exp_rk43_1st);
extern "C" define_e2v_slaveFunPtr(exp_rk43_1st);

/// The 2nd step of RK4 method
extern "C" define_e2v_FunPtr(exp_rk43_2nd);
extern "C" define_e2v_slaveFunPtr(exp_rk43_2nd);

/// The 3rd step of RK4 method
extern "C" define_e2v_FunPtr(exp_rk43_3rd);
extern "C" define_e2v_slaveFunPtr(exp_rk43_3rd);

/// The 4th step of RK4 method
extern "C" define_e2v_FunPtr(exp_rk43_4th);
extern "C" define_e2v_slaveFunPtr(exp_rk43_4th);

extern "C"
{

void exp_rk43_copy_host_(swInt* faceNum, swFloat* uc0, swFloat* uc)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 5, face_num, COPYOUT, uc0);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, uc);
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
 	opt.fun_slave = slave_exp_rk43_copy;

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

void exp_rk43_1st_host_(swInt* faceNum, swFloat* uc0, swFloat* uc,
	swFloat* rhs, swFloat* LHS_s, swFloat* is_int)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 5, face_num, COPYOUT, uc);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, uc0);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, rhs);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, LHS_s);
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
 	opt.fun_slave = slave_exp_rk43_1st;

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

void exp_rk43_2nd_host_(swInt* faceNum, swFloat* uc,
	swFloat* rhs, swFloat* LHS_s, swFloat* is_int)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 5, face_num, COPYINOUT, uc);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, rhs);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, LHS_s);
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
 	opt.fun_slave = slave_exp_rk43_2nd;

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

void exp_rk43_3rd_host_(swInt* faceNum, swFloat* uc, swFloat* uc0,
	swFloat* rhs, swFloat* LHS_s, swFloat* is_int, swFloat* R23,
	swFloat* R13, swFloat* R16)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 5, face_num, COPYINOUT, uc);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, uc0);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, rhs);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, LHS_s);
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
 	opt.fun_slave = slave_exp_rk43_3rd;

 	Arrays paraData;
 	swFloat para[3] = {*R23, *R13, *R16};
 	constructSingleArray(paraData, 1, 3, COPYIN, &para[0]);

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



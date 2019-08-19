#include <stdio.h>
#include <stdlib.h>
#include "swMacro.h"
#include "iterator.h"

#include "topology.hpp"
#include "rowSubsectionIterator.hpp"
// #include "wrappedInterface.h"

extern "C" define_e2v_FunPtr(gra_gg_div);
extern "C" define_e2v_slaveFunPtr(gra_gg_div);

extern "C" define_e2v_FunPtr(gra_gg_zero);
extern "C" define_e2v_slaveFunPtr(gra_gg_zero);

// 
extern "C"
{

 void gra_gg_host_(swInt* faceNum, swFloat* gra, swFloat* sec_vol)
 // void upw_conv_host_(swInt* faceNum)
 {
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 19, face_num, COPYINOUT, gra);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, sec_vol);
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
 	opt.fun_slave = slave_gra_gg_div;

 	Arrays paraData;
 	constructEmptyArray(paraData);

 	int cellWeights = 1;
 	int faceWeights = 20;
 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo;
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, ARRAY, false);
	iter->arrayIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
 }
//void upw_conv_host_(){}

void gra_gg_zero_host_(swInt* faceNum, swFloat* gra)
 // void upw_conv_host_(swInt* faceNum)
{
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 19, face_num, COPYINOUT, gra);
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
 	opt.fun_slave = slave_gra_gg_zero;

 	Arrays paraData;
 	constructEmptyArray(paraData);

 	int cellWeights = 1;
 	int faceWeights = 19;
 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo;
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, ARRAY, false);
	iter->arrayIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
}

}

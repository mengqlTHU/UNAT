#include <stdio.h>
#include <stdlib.h>
#include "swMacro.h"
#include "iterator.h"

#include "topology.hpp"
#include "rowSubsectionIterator.hpp"
// #include "wrappedInterface.h"

extern "C" define_e2v_FunPtr(uc2u);
extern "C" define_e2v_slaveFunPtr(uc2u);

// 
extern "C"
{

 void uc2u_host_(swInt* faceNum, swFloat* u, swFloat* uc, swFloat* t, 
 	swFloat* mu, swFloat* is_int, swFloat* rr, swFloat* gk1,
 	swFloat* t_suth, swFloat* is_tur_cal)
 // void upw_conv_host_(swInt* faceNum)
 {
 	swInt face_num = *faceNum;

 	Arrays frontEdgeData, backEdgeData, vertexData, selfConnData;
 	constructSingleArray(frontEdgeData, 5, face_num, COPYOUT, u);
 	addSingleArray(frontEdgeData, 5, face_num, COPYIN, uc);
 	addSingleArray(frontEdgeData, 1, face_num, COPYOUT, t);
 	addSingleArray(frontEdgeData, 1, face_num, COPYIN, is_int);
 	int cellWeights=1;
 	int faceWeights=1;
 	if(*is_tur_cal==1)
 	{
 		addSingleArray(frontEdgeData, 2, face_num, COPYINOUT, mu);
 	} else
 	{
 		addSingleArray(frontEdgeData, 1, face_num, COPYOUT, mu);
 	}
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
 	opt.fun_slave = slave_uc2u;

 	Arrays paraData;
 	swFloat para[3] = {*rr, *gk1, *t_suth};
 	constructSingleArray(paraData, 1, 3, COPYIN, &para[0]);

 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo;
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, ARRAY, false);
	iter->arrayIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
 }
//void upw_conv_host_(){}

}

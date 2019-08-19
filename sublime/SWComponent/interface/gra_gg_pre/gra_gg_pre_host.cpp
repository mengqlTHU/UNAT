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

extern "C" define_e2v_FunPtr(gra_gg_pre);
extern "C" define_e2v_slaveFunPtr(gra_gg_pre);
// extern UNAT::Iterator *iter;
// extern void constructRSSIterator(swInt* rowAddr, swInt* colAddr,
	// swInt* cellWeights, swInt* faceWeights, swInt faceNum);
// extern void constructMLBIterator(swInt* rowAddr, swInt* colAddr,
	// swInt* cellWeights, swInt* faceWeights, swInt faceNum);
// 
extern "C"
{

 void gra_gg_pre_host_(swInt* rowAddr, swInt* colAddr,
 	swInt* faceNum, swInt* cellNum,
 	swFloat* transform, swFloat* n_vg, swFloat* u, swFloat* gra, 
 	swFloat* sec_vol, swFloat* t, swFloat* is_int)
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
 	constructSingleArray(vertexData, 19, cell_num, COPYOUT, gra);
 	addSingleArray(vertexData, 5, cell_num, COPYIN, u);
 	addSingleArray(vertexData, 1, cell_num, COPYIN, sec_vol);
 	addSingleArray(vertexData, 1, cell_num, COPYIN, t);
 	addSingleArray(vertexData, 1, cell_num, COPYIN, is_int);
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
 	opt.fun_slave = slave_gra_gg_pre;

 	Arrays paraData;
 	constructEmptyArray(paraData);

 	int cellWeights = 27;
 	int faceWeights = 9;
 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo = UNAT::Topology::constructFromEdge(rowAddr, colAddr, face_num);
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, E2V, true);
	iter->edge2VertexIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
 }
//void upw_conv_host_(){}

}


// void upw_conv_host_(swInt* rowAddr, swInt* colAddr,
//  	swInt* cellWeights, swInt* faceWeights, swInt faceNum)
// {
// 	printf("%d\n",faceNum);
// 	// constructRSSIterator(rowAddr, colAddr, cellWeights, faceWeights, faceNum);
// }

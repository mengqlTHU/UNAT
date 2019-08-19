#include <stdio.h>
#include <stdlib.h>
#include "swMacro.h"
#include "iterator.h"

#include "topology.hpp"
#include "rowSubsectionIterator.hpp"
// #include "wrappedInterface.h"

extern "C" define_e2v_FunPtr(imp_lhs_imp);
extern "C" define_e2v_slaveFunPtr(imp_lhs_imp);

// 
extern "C"
{

 void imp_lhs_imp_host_(swInt* rowAddr, swInt* colAddr,
 	swInt* faceNum, swInt* cellNum, swFloat* transform, 
 	swFloat* n_vg, swFloat* sec_bct, swFloat* u, swFloat* mu,
 	swFloat* is_vis_cal, swFloat* is_tur_cal, swFloat* gk,
 	swFloat* gk1, swFloat* LHS_efix, swFloat* uc, swFloat* prl,
 	swFloat* prt, swFloat* rr, swFloat* cp)
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
 	constructSingleArray(vertexData, 5, cell_num, COPYIN, u);
 	addSingleArray(vertexData, 1, cell_num, COPYIN, sec_bct);
  	int cellWeights = 11;
 	int faceWeights = 9;
 	if(*is_tur_cal==1)
 	{
 		addSingleArray(vertexData, 2, cell_num, COPYIN, mu);
 		cellWeights += 2;
 	} else
 	{
 		addSingleArray(vertexData, 1, cell_num, COPYIN, mu);
 		cellWeights += 1;
 	}
 	addSingleArray(vertexData, 5, cell_num, COPYOUT, uc);
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
 	opt.fun_slave = slave_imp_lhs_imp;

 	Arrays paraData;
 	swFloat para[9] = {*is_vis_cal, *is_tur_cal, *gk, *gk1, *LHS_efix,
 		*prl, *prt, *rr, *cp};
 	constructSingleArray(paraData, 1, 9, COPYIN, &para[0]);


 	// constructRSSIterator(rowAddr, colAddr, &cellWeights, &faceWeights, face_num,
 		// &paraData, &opt, 1);
 	static UNAT::Topology *topo = UNAT::Topology::constructFromEdge(rowAddr, colAddr, face_num);
	static UNAT::Iterator *iter = new UNAT::RowSubsectionIterator(&opt, *topo, &cellWeights,
		&faceWeights, E2V, true);
	iter->edge2VertexIteration(&paraData, &opt, 1);
 	// constructMLBIterator(rowAddr, colAddr, cellWeights, faceWeights, *faceNum);
 }

}

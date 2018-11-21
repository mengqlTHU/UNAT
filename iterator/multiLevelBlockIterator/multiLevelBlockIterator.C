#include "multiLevelBlockIterator.H"
#include <stdlib.h>
#include <iostream>

#ifdef __cplusplus
extern "C"{
#endif
#include "BlockOrdering/BlockOrdering.h"
#include "edge2VertexIter_host.h"
#ifdef __cplusplus
}
#endif

#define MAX(A,B) A > B ? A : B

void MultiLevelBlockIterator::reorderEdgesFromEdge(swInt* startVertices,
			swInt* endVertices, swInt edgeNumber, swInt vertexNumber)
{
	swInt* vertexWeights;
	swInt* edgeWeights;
	swInt* postVertexOrder;
	swInt* postEdgeOrder;
	MLB_graph graph;
	swInt* blockNums;
	swInt  levels;

	postVertexOrder = (swInt*)malloc(sizeof(swInt)*vertexNumber);	
	postEdgeOrder = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	swInt totalSize
		// sendX, recvX, b
		= vertexNumber*3*sizeof(swFloat)
		// 对角块的owner，neighbour和非零元素
		+ edgeNumber*3/4*(2*sizeof(swInt)+sizeof(swFloat))
		// 非对角块的非零元素及owner（发送+接收）
		+ edgeNumber/4*2*(sizeof(swInt)+sizeof(swFloat))
		// 寄存器通信Pack（32Byte）
		+ edgeNumber/4/6*2*32
		// 发送接收计数数组和blockStart，cellStarts
		+ BLOCKNUM64K*BLOCKNUM64K*8*sizeof(swInt);

	this->_mshBlockNum = MAX(totalSize/64/64/1024+1,2);
	printf("The estimated mshBlockNum = %d\n",this->_mshBlockNum);

	vertexWeights = (swInt*)malloc(sizeof(swInt)*vertexNumber);
	edgeWeights = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	for(int i=0;i<edgeNumber;i++) edgeWeights[i] = 2;
	for(int i=0;i<vertexNumber;i++) vertexWeights[i] = 1;

	graph.owner       = startVertices;
	graph.neighbor    = endVertices;
	graph.cellWeights = vertexWeights;
	graph.edgeWeights = edgeWeights;
	graph.cellNum     = vertexNumber;
	graph.edgeNum     = edgeNumber;

	this->_cpeBlockNum = this->_mshBlockNum*BLOCKNUM64K;
	this->_mtxBlockNum = this->_cpeBlockNum*(this->_cpeBlockNum+1)/2;
	levels = 2;
	blockNums = (swInt*)malloc(sizeof(swInt)*levels);
	blockNums[0] = this->_mshBlockNum;
	blockNums[1] = BLOCKNUM64K;
	this->_blockStarts = (swInt*)malloc(sizeof(swInt)*this->_mtxBlockNum*4);
	this->_vertexStarts =
		(swInt*)malloc(sizeof(swInt)*(this->_cpeBlockNum+1));

	MLB_Multilevel_ordering(graph,levels,blockNums,this->_blockStarts,
				this->_vertexStarts,postVertexOrder,postEdgeOrder);

	this->_maxEdges = 0;
	this->_maxCells = 0;
	swInt edgeNum = 0;
	swInt cellNum = 0;
	swInt blockIdx;
	for(int i=0;i<this->_cpeBlockNum;i++)
	{
		for(int j=i+1;j<this->_cpeBlockNum;j++)
		{
			blockIdx = i*(1+2*this->_cpeBlockNum-i)/2+j-i;
			edgeNum += this->_blockStarts[4*blockIdx+3]
				- this->_blockStarts[4*blockIdx+2];
		}
		blockIdx = i*(1+2*this->_cpeBlockNum-i)/2;
		cellNum = this->_blockStarts[4*blockIdx+3]
			- this->_blockStarts[4*blockIdx+2];
		this->_maxCells
			= this->_maxCells > cellNum ? this->_maxCells : cellNum;
		this->_maxEdges
			= this->_maxEdges > edgeNum ? this->_maxEdges : edgeNum;
		edgeNum = 0;
	}
	this->_maxXNum = vertexNumber/this->_cpeBlockNum*1.2;

	this->_owner    = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	this->_neighbor = (swInt*)malloc(sizeof(swInt)*edgeNumber);

	for(int i=0;i<edgeNumber;i++)
	{
		this->_edgeMap.insert(pair<swInt,swInt>(i,postEdgeOrder[i]));
		if(postEdgeOrder[i]<0)
		{
			this->_owner[-postEdgeOrder[i]]
				= postVertexOrder[endVertices[i]];
			this->_neighbor[-postEdgeOrder[i]]
				= postVertexOrder[startVertices[i]];
		} else
		{
			this->_owner[postEdgeOrder[i]]
				= postVertexOrder[startVertices[i]];
			this->_neighbor[postEdgeOrder[i]]
				= postVertexOrder[endVertices[i]];
		}
	}
	for(int i=0;i<vertexNumber;i++)
	{
		this->_vertexMap.insert(pair<swInt,swInt>(i,postVertexOrder[i]));
	}

	cout<<"reorderEdgesFromEdge"<<endl;
}

void MultiLevelBlockIterator::reorderEdgesFromVertex(
			swInt* accuVertexEdgeNumbers, swInt* vertexNeighbours,
			swInt edgeNumber, swInt vertexNumber)
{
	cout<<"!!!---Not implemented in MLB"<<endl;
}

void MultiLevelBlockIterator::edge2VertexIteration(Arrays* edgeData,
			Arrays* vertexData,void (*operatorFunPointer)
			(MLBFunParameters *MLBFunParas))
{
	cout<<"operatorFunPointer"<<endl;

	MLBParameters MLBParas;
	MLBParas.blockStarts  = this->getBlockStarts();
	MLBParas.vertexStarts = this->getVertexStarts();
	MLBParas.owner        = this->getOwner();
	MLBParas.neighbor     = this->getNeighbor();
	MLBParas.cpeBlockNum  = this->getCpeBlockNum();
	MLBParas.mshBlockNum  = this->getMshBlockNum();
	MLBParas.mtxBlockNum  = this->getMtxBlockNum();
	MLBParas.maxXNum      = this->getMaxXNum();
	MLBParas.maxCells     = this->getMaxCells();
	MLBParas.maxEdges     = this->getMaxEdges();

//	edgeData->A3Ptr = this->getOwner();
//	edgeData->A4Ptr = this->getNeighbor();
//	operatorFunPointer(edgeData, vertexData);
	Arrays edgeData_slave, vertexData_slave;
	swInt edgeNum       = this->_topo->getEdgeNumber();
	swInt vertexNum     = this->_topo->getVertexNumber();
	swInt cpeBlockNum   = this->getCpeBlockNum();
	swInt mshBlockNum   = this->getMshBlockNum();
	swInt mtxBlockNum   = this->getMtxBlockNum();
	swInt* blockStarts  = this->getBlockStarts();
	swInt* vertexStarts = this->getVertexStarts();
	swInt* owner        = this->getOwner();
	swInt* neighbor     = this->getNeighbor();

//	for(int i=0;i<vertexNum;i++)
//	{
//		vertexData->A1Ptr[i] = vertexData->A2Ptr[i]*vertexData->A3Ptr[i];
//	}

	edge2VertexIteration_host(edgeData, vertexData,
				operatorFunPointer, &MLBParas);

//	for(int i=0;i<mshBlockNum;i++)
//	{
//		for(int j=0;j<BLOCKNUM64K;j++)
//		{
//			int row = i*BLOCKNUM64K+j;
//			for(int col=row;col<cpeBlockNum;col++)
//			{
//				int blockIdx = row*(1+2*cpeBlockNum-row)/2+col-row;
//				int startIdx = blockStarts[4*blockIdx+2];
//				int endIdx   = blockStarts[4*blockIdx+3];
//				for(int k=startIdx;k<endIdx;k++)
//				{
//				edgeData_slave.A2Ptr = &edgeData->A2Ptr[k];
//				vertexData_slave.A1Ptr = &vertexData->A1Ptr[owner[k]];
//				vertexData_slave.A2Ptr = &vertexData->A2Ptr[neighbor[k]];
//				operatorFunPointer(&edgeData_slave, &vertexData_slave,
//							&owner[k], &neighbor[k]);
//				edgeData_slave.A2Ptr = &edgeData->A1Ptr[k];
//				vertexData_slave.A1Ptr = &vertexData->A1Ptr[neighbor[k]];
//				vertexData_slave.A2Ptr = &vertexData->A2Ptr[owner[k]];
//				operatorFunPointer(&edgeData_slave, &vertexData_slave,
//							&owner[k], &neighbor[k]);
////					if(owner[k]==30527)
////					{
////						cout<<"row: "<<k<<","
////							<<vertexData->A1Ptr[owner[k]]<<","
////							<<edgeData->A2Ptr[k]<<","
////							<<vertexData->A2Ptr[neighbor[k]]<<endl;
////					}
////					vertexData->A1Ptr[owner[k]]
////						+= edgeData->A2Ptr[k]
////						 * vertexData->A2Ptr[neighbor[k]];
////					if(neighbor[k]==30534)
////					{
////						cout<<"col: "<<k<<","
////							<<vertexData->A1Ptr[neighbor[k]]<<","
////							<<edgeData->A1Ptr[k]<<","
////							<<vertexData->A2Ptr[owner[k]]<<endl;;
////					}
////					vertexData->A1Ptr[neighbor[k]]
////						+= edgeData->A1Ptr[k]
////						 * vertexData->A2Ptr[owner[k]];
//				}
//			}
//		}
//	}

	map<swInt, swInt>::iterator iter;
	vertexData->A4Ptr=(swFloat*)malloc
		(this->_topo->getVertexNumber()*sizeof(swFloat));
	for(iter = this->_vertexMap.begin();iter!=this->_vertexMap.end();iter++)
	{
		vertexData->A4Ptr[iter->first]=vertexData->A1Ptr[iter->second];
	}
}

void MultiLevelBlockIterator::vertex2EdgeIteration(Arrays* neighbourData,
			Arrays* vertexData,void (*operatorFunPointer)
			(MLBFunParameters *MLBFunParas))
{
	cout<<"!!!---Not implemented in MLB"<<endl;
//	operatorFunPointer(edgeData, vertexData, this->_topo);
}

void MultiLevelBlockIterator::reorderEdgeData(Arrays* edgeData)
{
	swFloat* lower = edgeData->A1Ptr;
	swFloat* upper = edgeData->A2Ptr;
	cout<<"reorderEdgeData"<<endl;
	map<swInt, swInt>::iterator iter;
	swFloat* tmpL=(swFloat*)malloc
		(this->_topo->getEdgeNumber()*sizeof(swFloat));
	swFloat* tmpU=(swFloat*)malloc
		(this->_topo->getEdgeNumber()*sizeof(swFloat));
	for(iter = this->_edgeMap.begin();iter!=_edgeMap.end();iter++)
	{
		if(iter->second<0)
		{
			tmpL[-iter->second] = upper[iter->first];
			tmpU[-iter->second] = lower[iter->first];
		} else
		{
			tmpL[iter->second] = lower[iter->first];
			tmpU[iter->second] = upper[iter->first];
		}
	}
	for(int i=0;i<this->_topo->getEdgeNumber();i++)
	{
		edgeData->A1Ptr[i] = tmpL[i];
		edgeData->A2Ptr[i] = tmpU[i];
	}
	free(tmpL);
	free(tmpU);
}

void MultiLevelBlockIterator::reorderVertexArray(swFloat* array)
{
	swFloat *tmp = (swFloat*)malloc
		(this->_topo->getVertexNumber()*(sizeof(swFloat)));
	map<swInt, swInt>::iterator iter;
	for(iter = this->_vertexMap.begin();iter!=this->_vertexMap.end();iter++)
	{
//		if(iter->first==0) cout<<iter->second<<endl;
		tmp[iter->second] = array[iter->first];
	}
	for(int i=0;i<this->_topo->getVertexNumber();i++)
	{
		array[i] = tmp [i];
	}
	free(tmp);
}

void MultiLevelBlockIterator::reorderVertexData(Arrays* vertexData)
{
	if(vertexData->A1Ptr!=NULL) reorderVertexArray(vertexData->A1Ptr);
	if(vertexData->A2Ptr!=NULL) reorderVertexArray(vertexData->A2Ptr);
	if(vertexData->A3Ptr!=NULL) reorderVertexArray(vertexData->A3Ptr);
	if(vertexData->A4Ptr!=NULL) reorderVertexArray(vertexData->A4Ptr);
}


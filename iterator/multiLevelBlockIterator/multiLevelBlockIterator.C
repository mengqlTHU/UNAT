#include "multiLevelBlockIterator.H"
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#ifdef __cplusplus
extern "C"{
#endif
#include "BlockOrdering/BlockOrdering.h"
#include "edge2VertexIter_host.h"
#include "vertex2EdgeIter_host.h"
#ifdef __cplusplus
}
#endif

#define MAX(A,B) A > B ? A : B

MultiLevelBlockIterator::MultiLevelBlockIterator(Topology &topo)
		: Iterator(topo)
{
//	LOG("MultiLevelBlockIterator");
    printf("start MultiLevelBlocking reorder...\n");
	MLBReorder(topo,2);
	printf("The domain is divided into %d*%d blocks.\n",
				this->_mshBlockNum,this->_mshBlockNum);

	swInt vertexNumber = this->getTopology()->getVertexNumber();
	swInt edgeNumber   = this->getTopology()->getEdgeNumber();

	for(int i=0;i<edgeNumber;i++)
	{
		this->getEdgeMap()
			.insert(pair<swInt,swInt>(i,this->_postEdgeOrder[i]));
	}
	for(int i=0;i<vertexNumber;i++)
	{
		this->getVertexMap()
			.insert(pair<swInt,swInt>(i,this->_postVertexOrder[i]));
	}
}

void MultiLevelBlockIterator::MLBReorder(Topology &topo, swInt ref)
{
	swInt* vertexWeights;
	swInt* edgeWeights;
	MLB_graph graph;
	swInt* blockNums;
	swInt  levels;
	swInt vertexNumber = this->getTopology()->getVertexNumber();
	swInt edgeNumber   = this->getTopology()->getEdgeNumber();

	this->_postVertexOrder = (swInt*)malloc(sizeof(swInt)*vertexNumber);	
	this->_postEdgeOrder = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	swInt* postVertexOrder = this->_postVertexOrder;;
	swInt* postEdgeOrder = this->_postEdgeOrder;;
	swInt totalSize
		// sendX, recvX, b
		= vertexNumber*2*sizeof(swFloat)
		// 对角块的owner，neighbour和非零元素
		+ edgeNumber*2*2*(sizeof(swInt)+sizeof(swFloat))
		// 发送接收计数数组和blockStart，cellStarts
		+ BLOCKNUM64K*BLOCKNUM64K*5*sizeof(swInt);

	this->_mshBlockNum = MAX(totalSize/64/64/1024+1,25);
	printf("The estimated mshBlockNum = %d\n",this->_mshBlockNum);

	vertexWeights = (swInt*)malloc(sizeof(swInt)*vertexNumber);
	edgeWeights = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	for(int i=0;i<edgeNumber;i++) edgeWeights[i] = 2;
	for(int i=0;i<vertexNumber;i++) vertexWeights[i] = 1;

	graph.owner       = this->getTopology()->getStartVertices();
	graph.neighbor    = this->getTopology()->getEndVertices();
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
	this->_maxXNum  = 0;
	swInt xNum    = 0;
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
		xNum = this->_vertexStarts[i+1]-this->_vertexStarts[i];
		this->_maxXNum
			= this->_maxXNum > xNum ? this->_maxXNum : xNum;	
		this->_maxCells
			= this->_maxCells > cellNum ? this->_maxCells : cellNum;
		this->_maxEdges
			= this->_maxEdges > edgeNum ? this->_maxEdges : edgeNum;
		edgeNum = 0;
	}
	this->_maxEdgesUnsymm = (this->_maxCells*2+this->_maxEdges)*1.01;
	swInt totalLength
		= this->_maxXNum*2*sizeof(swFloat)
		+ this->_maxEdgesUnsymm*2*sizeof(swFloat)
		+ this->_maxEdgesUnsymm*2*sizeof(swInt)
		+ BLOCKNUM64K*5*sizeof(swInt)
		+ this->_maxEdges/6*2*32;
//	if(totalLength > 64*1024*0.97) MLBReorder(topo,this->_mshBlockNum+1);
}

void MultiLevelBlockIterator::reorderEdgesFromEdge(swInt* startVertices,
			swInt* endVertices, swInt edgeNumber, swInt vertexNumber)
{
	if(edgeNumber!=this->getTopology()->getEdgeNumber())
	{
		LOG("The edge arrays do not match the topology!");
	}

	map<swInt, swInt>::iterator iter;
	this->_owner   = (swInt*)malloc(sizeof(swInt)*edgeNumber);
	this->_neighbor = (swInt*)malloc(sizeof(swInt)*edgeNumber);

	for(iter = this->getEdgeMap().begin();
				iter!=this->getEdgeMap().end();iter++)
	{
		if(iter->second<0)
		{
			this->_owner[-iter->second]
				= this->getVertexMap()[endVertices[iter->first]];
			this->_neighbor[-iter->second]
				= this->getVertexMap()[startVertices[iter->first]];
		} else
		{
			this->_owner[iter->second]
				= this->getVertexMap()[startVertices[iter->first]];
			this->_neighbor[iter->second]
				= this->getVertexMap()[endVertices[iter->first]];
		}
	}

	for(int i=0;i<edgeNumber;i++)
	{
		startVertices[i] = this->_owner[i];
		endVertices[i]   = this->_neighbor[i];
	}
	
//	LOG("reorderEdgesFromEdge");
}

void MultiLevelBlockIterator::reorderEdgesFromVertex(
			swInt* firstEdgeVertices, swInt* vertexNeighbours,
			swInt edgeNumber, swInt vertexNumber)
{
	this->_firstEdgeVertices
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	this->_vertexNeighbours
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	for(int i=0;i<edgeNumber*2;i++) 
	{
		this->_firstEdgeVertices[i] = firstEdgeVertices[i];
		this->_vertexNeighbours[i]  = vertexNeighbours[i];
	}

	swInt cpeBlockNum   = this->getCpeBlockNum();
	swInt mshBlockNum   = this->getMshBlockNum();
	swInt mtxBlockNum   = this->getMtxBlockNum();
	swInt* blockStarts  = this->getBlockStarts();
	swInt* vertexStarts = this->getVertexStarts();
	swInt* owner        = this->getTopology()->getStartVertices();
	swInt* neighbor     = this->getTopology()->getEndVertices();

	swInt* vertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	swInt* accuVertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	swInt* tmpVN
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	swInt* tmpFEV
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	this->_blockStartsUnsymm
		= (swInt*)malloc(sizeof(swInt)*cpeBlockNum*cpeBlockNum*4);
	swInt* blockStartsUnsymm = this->_blockStartsUnsymm;

	int blockIdx, startBlockIdx, endBlockIdx, row1, col1, blockIdxUnsymm;
	int accuEdges= 0;
	int startIdx = 0;
	int endIdx   = 0;
	for(int i=0;i<mshBlockNum;i++)
	{
		for(int j=0;j<BLOCKNUM64K;j++)
		{
			int row = i*BLOCKNUM64K+j;
			int cellLen = vertexStarts[row+1]-vertexStarts[row];
			assert(cellLen<vertexNumber/cpeBlockNum*2);
			for(int col=0;col<cpeBlockNum;col++)
			{
				blockIdxUnsymm = row*cpeBlockNum+col;
				blockStartsUnsymm[4*blockIdxUnsymm  ] = row;	
				blockStartsUnsymm[4*blockIdxUnsymm+1] = col;	
				for(int k=0;k<vertexNumber/cpeBlockNum*2;k++)
				{
					vertexEdgeNumbers[k]=0;
					accuVertexEdgeNumbers[k]=0;
				}
				if(col==row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
		    			vertexEdgeNumbers[owner[k]-vertexStarts[row]]++;
		    			vertexEdgeNumbers[neighbor[k]-vertexStarts[row]]++;
		    		}
					for(int k=1;k<cellLen;k++)
					{
						assert(vertexEdgeNumbers[k-1]>=0);
						accuVertexEdgeNumbers[k+1]
							= accuVertexEdgeNumbers[k]
							+ vertexEdgeNumbers[k-1];
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = neighbor[k]+1-vertexStarts[row];
						col1 = owner[k];
						tmpVN[accuVertexEdgeNumbers[row1]+startIdx]=col1;
						tmpFEV[accuVertexEdgeNumbers[row1]+startIdx]
							=neighbor[k];
						accuVertexEdgeNumbers[row1]++;
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = owner[k]+1-vertexStarts[row];
						col1 = neighbor[k];
						tmpVN[accuVertexEdgeNumbers[row1]+startIdx]=col1;
						tmpFEV[accuVertexEdgeNumbers[row1]+startIdx]
							=owner[k];
						accuVertexEdgeNumbers[row1]++;
					}
					blockStartsUnsymm[4*blockIdxUnsymm+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx)*2;
					blockStartsUnsymm[4*blockIdxUnsymm+3] = startIdx;	
				} else if(col < row)
				{
     				blockIdx = col*(1+2*this->_cpeBlockNum-col)/2+row-col;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
						tmpVN[k-startBlockIdx+startIdx] = owner[k];
						tmpFEV[k-startBlockIdx+startIdx] = neighbor[k];
		    		}
					blockStartsUnsymm[4*blockIdxUnsymm+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx);
					blockStartsUnsymm[4*blockIdxUnsymm+3] = startIdx;	
				} else if(col > row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2+col-row;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
						tmpVN[k-startBlockIdx+startIdx] = neighbor[k];
						tmpFEV[k-startBlockIdx+startIdx] = owner[k];
		    		}
					blockStartsUnsymm[4*blockIdxUnsymm+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx);
					blockStartsUnsymm[4*blockIdxUnsymm+3] = startIdx;	
				}
			}
		}
	}
	for(int i=0;i<edgeNumber*2;i++)
	{
		firstEdgeVertices[i] = tmpFEV[i];
		vertexNeighbours[i] = tmpVN[i];
	}
	swInt maxEdgesUnsymm = 0;
	swInt edgeNum = 0;
	for(int i=0;i<cpeBlockNum;i++)
	{
		edgeNum = blockStartsUnsymm[4*(i+1)*(cpeBlockNum-1)+3]
			- blockStartsUnsymm[4*i*cpeBlockNum+2];
		maxEdgesUnsymm = maxEdgesUnsymm>edgeNum?maxEdgesUnsymm:edgeNum;
	}
	assert(this->_maxEdgesUnsymm>maxEdgesUnsymm);
	free(tmpFEV);
	free(tmpVN);
	free(vertexEdgeNumbers);
	free(accuVertexEdgeNumbers);
//	assert(startIdx==edgeNumber*2);
//	for(int i=0;i<edgeNumber*2;i++)
//	{
//		if(vertexNeighbours[i]==-1 || firstEdgeVertices[i]==-1)
//		{
//			printf("%d,%d,%d\n",i,vertexNeighbours[i],firstEdgeVertices[i]);
//			exit(1);
//		}
//	}
//
//	for(int i=0;i<mshBlockNum;i++)
//	{
//		for(int j=0;j<BLOCKNUM64K;j++)
//		{
//			int row = i*BLOCKNUM64K+j;
//			for(int col=0;col<cpeBlockNum;col++)
//			{
//				blockIdxV = row*cpeBlockNum+col;
//				printf("(%4d,%4d),",blockStartsV[4*blockIdxV+2],blockStartsV[4*blockIdxV+3]);
//			}
//		}
//	}
//printf("row:%d, col:%d, edgeNums:%d, vertexNums:%d, startIdx:%d\n",row,col,endBlockIdx-startBlockIdx,cellLen,startIdx);
//for(int m=0;m<140;m++)
//{
//	for(int n=0;n<15;n++)
//	{
//		printf("%4d ",owner[m*15+n]);
//	}
//	printf("\n");
//}
//printf("********************************************\n");
//	LOG("reorderEdgesFromVertex");
}

void MultiLevelBlockIterator::edge2VertexIteration(Arrays* edgeData,
			Arrays* vertexData,void(*operatorFunPointer_host)
			(MLBFunParameters *MLBFunParas), void(*operatorFunPointer_slave)
			(MLBFunParameters *MLBFunParas))
{
//	LOG("operatorFunPointer");

	MLBParameters MLBParas;
	MLBParas.blockStarts  = this->getBlockStarts();
	MLBParas.vertexStarts = this->getVertexStarts();
	MLBParas.owner        = this->getTopology()->getStartVertices();
	MLBParas.neighbor     = this->getTopology()->getEndVertices();
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
	swInt edgeNum       = this->getTopology()->getEdgeNumber();
	swInt vertexNum     = this->getTopology()->getVertexNumber();
	swInt cpeBlockNum   = this->getCpeBlockNum();
	swInt mshBlockNum   = this->getMshBlockNum();
	swInt mtxBlockNum   = this->getMtxBlockNum();
	swInt* blockStarts  = this->getBlockStarts();
	swInt* vertexStarts = this->getVertexStarts();
	swInt* owner        = this->getTopology()->getStartVertices();
	swInt* neighbor     = this->getTopology()->getEndVertices();
	MLBParas.operatorFunPointer_host  = operatorFunPointer_host;
	MLBParas.operatorFunPointer_slave  = operatorFunPointer_slave;
	swInt dimension     = edgeData->dimension;

//	for(int i=0;i<vertexNum;i++)
//	{
//		vertexData->A1Ptr[i] = vertexData->A2Ptr[i]*vertexData->A3Ptr[i];
//	}

	edge2VertexIteration_host(edgeData, vertexData, &MLBParas);

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
		(this->getTopology()->getVertexNumber()*sizeof(swFloat)*dimension);
	for(int j=0;j<dimension;j++)
	{
		for(iter = this->getVertexMap().begin();
					iter!=this->getVertexMap().end();iter++)
		{
			vertexData->A4Ptr[iter->first+j*vertexNum]
				= vertexData->A1Ptr[iter->second+j*vertexNum];
		}
	}
}

void MultiLevelBlockIterator::vertex2EdgeIteration(Arrays* neighbourData,
			Arrays* vertexData,void(*operatorFunPointer_host)
			(MLBFunParameters *MLBFunParas), void(*operatorFunPointer_slave)
			(MLBFunParameters *MLBFunParas))
{
	MLBParameters MLBParas;
	MLBParas.blockStarts       = this->getBlockStarts();
	MLBParas.blockStartsUnsymm = this->getBlockStartsUnsymm();
	MLBParas.vertexStarts      = this->getVertexStarts();
	MLBParas.firstEdgeVertices
		= this->getTopology()->getFirstEdgeVertices();
	MLBParas.vertexNeighbor
		= this->getTopology()->getVertexNeighbours();
	MLBParas.owner        = this->getTopology()->getStartVertices();
	MLBParas.neighbor     = this->getTopology()->getEndVertices();
	MLBParas.cpeBlockNum  = this->getCpeBlockNum();
	MLBParas.mshBlockNum  = this->getMshBlockNum();
	MLBParas.mtxBlockNum  = this->getMtxBlockNum();
	MLBParas.maxXNum      = this->getMaxXNum();
	MLBParas.maxCells     = this->getMaxCells();
	MLBParas.maxEdges     = this->getMaxEdgesUnsymm();
	MLBParas.operatorFunPointer_host  = operatorFunPointer_host;
	MLBParas.operatorFunPointer_slave = operatorFunPointer_slave;

	vertex2EdgeIteration_host(neighbourData, vertexData, &MLBParas);

	map<swInt, swInt>::iterator iter;
	vertexData->A4Ptr=(swFloat*)malloc
		(this->getTopology()->getVertexNumber()*sizeof(swFloat));
	for(iter = this->getVertexMap().begin();
				iter!=this->getVertexMap().end();iter++)
	{
//		if(iter->first==0) printf("%d\n",iter->second);
		vertexData->A4Ptr[iter->first]=vertexData->A1Ptr[iter->second];
	}

//	LOG("vertex2EdgeIteration");
//	operatorFunPointer(edgeData, vertexData, this->_topo);
}

void MultiLevelBlockIterator::reorderEdgeDataUnsymm(Arrays* edgeData)
{
	if(edgeData->A1Ptr!=NULL) reorderEdgeArrayUnsymm(edgeData->A1Ptr);
	if(edgeData->A2Ptr!=NULL) reorderEdgeArrayUnsymm(edgeData->A2Ptr);
	if(edgeData->A3Ptr!=NULL) reorderEdgeArrayUnsymm(edgeData->A3Ptr);
	if(edgeData->A4Ptr!=NULL) reorderEdgeArrayUnsymm(edgeData->A4Ptr);
}

void MultiLevelBlockIterator::reorderEdgeArrayUnsymm(swFloat* array)
{
	swInt  edgeNumber    = this->getTopology()->getEdgeNumber();
	swInt  vertexNumber  = this->getTopology()->getVertexNumber();
	swInt  cpeBlockNum   = this->_cpeBlockNum;
	swInt  mshBlockNum   = this->_mshBlockNum;
	swInt* blockStarts  = this->getBlockStarts();
	swInt* vertexStarts = this->getVertexStarts();
	swInt* owner        = this->getTopology()->getStartVertices();
	swInt* neighbor     = this->getTopology()->getEndVertices();

	swFloat* tmp = (swFloat*)malloc(sizeof(swFloat)*edgeNumber*2);
	swInt* firstEdgeVertices = this->_firstEdgeVertices;
	swInt* vertexNeighbours  = this->_vertexNeighbours;
	swInt* accuVertexEdgeNumbers
		= this->getTopology()->getAccuVertexEdgeNumbers();
	swInt* accuStartVertexNumbers
		= this->getTopology()->getAccuStartVertexNumbers();
	swInt row,col,upIdx=0;
	swInt*  lowIdx = (swInt*)malloc(sizeof(swInt)*vertexNumber);
	swFloat* lower = (swFloat*)malloc(sizeof(swFloat)*edgeNumber);
	swFloat* upper = (swFloat*)malloc(sizeof(swFloat)*edgeNumber);
	for(int i=0;i<vertexNumber;i++) {lowIdx[i]=0;}
	for(int i=0;i<edgeNumber;i++) {lower[i]=-1;upper[i]=-1;}

	// COO->LDU
	for(int i=0;i<vertexNumber;i++)
	{
		for(int j=accuVertexEdgeNumbers[i];j<accuVertexEdgeNumbers[i+1];j++)
		{
			if(vertexNeighbours[j]>i)
			{
				upper[upIdx] = array[j];
				upIdx++;
			}else if(vertexNeighbours[j]<i)
			{
				assert(accuStartVertexNumbers[vertexNeighbours[j]]+lowIdx[vertexNeighbours[j]]<accuStartVertexNumbers[vertexNeighbours[j]+1]);
				lower[accuStartVertexNumbers[vertexNeighbours[j]]
					+lowIdx[vertexNeighbours[j]]] = array[j];
				lowIdx[vertexNeighbours[j]]++;
			}
		}
	}
//	for(int i=0;i<edgeNumber;i++)
//	{
//		if(upper[i]==-1 || lower[i]==-1)
//		{
//			printf("%d,%f,%f\n",i,lower[i],upper[i]);
//		}
//	}

	// LDU->reordered LDU
	Arrays edgeData = {lower, upper, NULL, NULL, edgeNumber};
	reorderEdgeData(&edgeData);

	// reordered LDU -> reordered COO data
	swInt* vertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	swInt* accuVertexEdgeNums
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);

	int blockIdx, startBlockIdx, endBlockIdx, row1, blockIdxV;
	swFloat value;
	int accuEdges= 0;
	int startIdx = 0;
	int endIdx   = 0;
	for(int i=0;i<mshBlockNum;i++)
	{
		for(int j=0;j<BLOCKNUM64K;j++)
		{
			int row = i*BLOCKNUM64K+j;
			int cellLen = vertexStarts[row+1]-vertexStarts[row];
			assert(cellLen<vertexNumber/cpeBlockNum*2);
			for(int col=0;col<cpeBlockNum;col++)
			{
				for(int k=0;k<vertexNumber/cpeBlockNum*2;k++)
				{
					vertexEdgeNumbers[k]=0;
					accuVertexEdgeNums[k]=0;
				}
				if(col==row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
		    			vertexEdgeNumbers[owner[k]-vertexStarts[row]]++;
		    			vertexEdgeNumbers[neighbor[k]-vertexStarts[row]]++;
		    		}
					for(int k=1;k<cellLen;k++)
					{
						accuVertexEdgeNums[k+1]
							= accuVertexEdgeNums[k]
							+ vertexEdgeNumbers[k-1];
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = neighbor[k]+1-vertexStarts[row];
//						if(accuVertexEdgeNums[row1]+startIdx==178894) printf("%d,%f\n",k,lower[k]);
						array[accuVertexEdgeNums[row1]+startIdx]=lower[k];
						accuVertexEdgeNums[row1]++;
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = owner[k]+1-vertexStarts[row];
//						if(accuVertexEdgeNums[row1]+startIdx==178894) printf("%d,%f\n",k,upper[k]);
						array[accuVertexEdgeNums[row1]+startIdx]=upper[k];
						accuVertexEdgeNums[row1]++;
					}
					startIdx += (endBlockIdx-startBlockIdx)*2;
				} else if(col < row)
				{
     				blockIdx = col*(1+2*this->_cpeBlockNum-col)/2+row-col;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
//						if(k-startBlockIdx+startIdx==178894) printf("%d,%f\n",k,lower[k]);
						array[k-startBlockIdx+startIdx] = lower[k];
		    		}
					startIdx += (endBlockIdx-startBlockIdx);
				} else if(col > row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2+col-row;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
//						if(k-startBlockIdx+startIdx==178894) printf("%d,%f\n",k,upper[k]);
						array[k-startBlockIdx+startIdx] = upper[k];
		    		}
					startIdx += (endBlockIdx-startBlockIdx);
				}
			}
		}
	}

	free(tmp);
	free(lower);
	free(upper);
	free(lowIdx);
	free(vertexEdgeNumbers);
	free(accuVertexEdgeNums);
}


void MultiLevelBlockIterator::reorderEdgeData(Arrays* edgeData)
{
	swFloat* lower = edgeData->A1Ptr;
	swFloat* upper = edgeData->A2Ptr;
//	LOG("reorderEdgeData");
	map<swInt, swInt>::iterator iter;
	swInt dimension = edgeData->dimension;
	swInt num       = edgeData->num;
	swFloat* tmpL
		= (swFloat*)malloc(edgeData->num*sizeof(swFloat)*dimension);
	swFloat* tmpU
		= (swFloat*)malloc(edgeData->num*sizeof(swFloat)*dimension);
	for(int j=0;j<dimension;j++)
	{
		for(iter = this->getEdgeMap().begin();
					iter!=this->getEdgeMap().end();iter++)
		{
			if(iter->second<0)
			{
				tmpL[-iter->second+j*num] = upper[iter->first+j*num];
				tmpU[-iter->second+j*num] = lower[iter->first+j*num];
			} else
			{
				tmpL[iter->second+j*num] = lower[iter->first+j*num];
				tmpU[iter->second+j*num] = upper[iter->first+j*num];
			}
		}
		for(int i=0;i<this->getTopology()->getEdgeNumber();i++)
		{
			edgeData->A1Ptr[i+j*num] = tmpL[i+j*num];
			edgeData->A2Ptr[i+j*num] = tmpU[i+j*num];
		}
	}
	free(tmpL);
	free(tmpU);
}

void MultiLevelBlockIterator::reorderVertexArray(swFloat* array)
{
	swFloat *tmp = (swFloat*)malloc
		(this->getTopology()->getVertexNumber()*(sizeof(swFloat)));
	map<swInt, swInt>::iterator iter;
	for(iter = this->getVertexMap().begin();
				iter!=this->getVertexMap().end();iter++)
	{
//		if(iter->first==0) cout<<iter->second<<endl;
		tmp[iter->second] = array[iter->first];
	}
	for(int i=0;i<this->getTopology()->getVertexNumber();i++)
	{
		array[i] = tmp [i];
	}
	free(tmp);
}

void MultiLevelBlockIterator::reorderVertexData(Arrays* vertexData)
{
	swInt dimension = vertexData->dimension;
	for(int j=0;j<dimension;j++)
	{
		if(vertexData->A1Ptr!=NULL) 
		  reorderVertexArray(&vertexData->A1Ptr[j*vertexData->num]);
		if(vertexData->A2Ptr!=NULL) 
		  reorderVertexArray(&vertexData->A2Ptr[j*vertexData->num]);
		if(vertexData->A3Ptr!=NULL)
		  reorderVertexArray(&vertexData->A3Ptr[j*vertexData->num]);
		if(vertexData->A4Ptr!=NULL)
		  reorderVertexArray(&vertexData->A4Ptr[j*vertexData->num]);
	}
}


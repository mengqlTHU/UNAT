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
	LOG("MultiLevelBlockIterator");
	swInt* vertexWeights;
	swInt* edgeWeights;
	swInt* postVertexOrder;
	swInt* postEdgeOrder;
	MLB_graph graph;
	swInt* blockNums;
	swInt  levels;
	swInt vertexNumber = this->getTopology()->getVertexNumber();
	swInt edgeNumber   = this->getTopology()->getEdgeNumber();

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

	for(int i=0;i<edgeNumber;i++)
	{
		this->getEdgeMap().insert(pair<swInt,swInt>(i,postEdgeOrder[i]));
	}
	for(int i=0;i<vertexNumber;i++)
	{
		this->getVertexMap()
			.insert(pair<swInt,swInt>(i,postVertexOrder[i]));
	}

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
	
	LOG("reorderEdgesFromEdge");
}

void MultiLevelBlockIterator::reorderEdgesFromVertex(
			swInt* firstEdgeVertices, swInt* vertexNeighbours,
			swInt edgeNumber, swInt vertexNumber)
{
	swInt cpeBlockNum   = this->getCpeBlockNum();
	swInt mshBlockNum   = this->getMshBlockNum();
	swInt mtxBlockNum   = this->getMtxBlockNum();
	swInt* blockStarts  = this->getBlockStarts();
	swInt* blockStartsV = this->getBlockStartsV();
	swInt* vertexStarts = this->getVertexStarts();
	swInt* owner        = this->getTopology()->getStartVertices();
	swInt* neighbor     = this->getTopology()->getEndVertices();

	swInt* vertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	swInt* accuVertexEdgeNumbers
		= (swInt*)malloc(sizeof(swInt)*vertexNumber/cpeBlockNum*2);
	vertexNeighbours
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	firstEdgeVertices
		= (swInt*)malloc(sizeof(swInt)*edgeNumber*2);
	blockStartsV
		= (swInt*)malloc(sizeof(swInt)*cpeBlockNum*cpeBlockNum*4);

//	for(int i=0;i<edgeNumber*2;i++) {vertexNeighbours[i]=-1;}
//	for(int i=0;i<edgeNumber*2;i++) {firstEdgeVertices[i]=-1;}

	int blockIdx, startBlockIdx, endBlockIdx, row1, col1, blockIdxV;
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
				blockIdxV = row*cpeBlockNum+col;
				blockStartsV[4*blockIdxV  ] = row;	
				blockStartsV[4*blockIdxV+1] = col;	
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
						vertexNeighbours
							[accuVertexEdgeNumbers[row1]+startIdx]=col1;
						firstEdgeVertices
							[accuVertexEdgeNumbers[row1]+startIdx]
							=neighbor[k];
						accuVertexEdgeNumbers[row1]++;
					}
					for(int k=startBlockIdx;k<endBlockIdx;k++)
					{
						row1 = owner[k]+1-vertexStarts[row];
						col1 = neighbor[k];
						vertexNeighbours
							[accuVertexEdgeNumbers[row1]+startIdx]=col1;
						firstEdgeVertices
							[accuVertexEdgeNumbers[row1]+startIdx]
							=owner[k];
						accuVertexEdgeNumbers[row1]++;
					}
					blockStartsV[4*blockIdxV+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx)*2;
					blockStartsV[4*blockIdxV+3] = startIdx;	
				} else if(col < row)
				{
     				blockIdx = col*(1+2*this->_cpeBlockNum-col)/2+row-col;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
						vertexNeighbours[k-startBlockIdx+startIdx]
							= owner[k];
						firstEdgeVertices[k-startBlockIdx+startIdx]
							= neighbor[k];
		    		}
					blockStartsV[4*blockIdxV+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx);
					blockStartsV[4*blockIdxV+3] = startIdx;	
				} else if(col > row)
				{
     				blockIdx = row*(1+2*this->_cpeBlockNum-row)/2+col-row;
     				startBlockIdx = this->_blockStarts[4*blockIdx+2];
     				endBlockIdx   = this->_blockStarts[4*blockIdx+3];
		    		for(int k=startBlockIdx;k<endBlockIdx;k++)
		    		{
						vertexNeighbours[k-startBlockIdx+startIdx]
							= neighbor[k];
						firstEdgeVertices[k-startBlockIdx+startIdx]
							= owner[k];
		    		}
					blockStartsV[4*blockIdxV+2] = startIdx;	
					startIdx += (endBlockIdx-startBlockIdx);
					blockStartsV[4*blockIdxV+3] = startIdx;	
				}
			}
		}
	}
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
//for(int m=0;m<40;m++)
//{
//	for(int n=0;n<15;n++)
//	{
//		printf("%4d ",accuVertexEdgeNumbers[m*15+n]);
//	}
//	printf("\n");
//}
//printf("********************************************\n");
	LOG("reorderEdgesFromVertex");
}

void MultiLevelBlockIterator::edge2VertexIteration(Arrays* edgeData,
			Arrays* vertexData,void(*operatorFunPointer_host)
			(MLBFunParameters *MLBFunParas), void(*operatorFunPointer_slave)
			(MLBFunParameters *MLBFunParas))
{
	LOG("operatorFunPointer");

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

//	for(int i=0;i<vertexNum;i++)
//	{
//		vertexData->A1Ptr[i] = vertexData->A2Ptr[i]*vertexData->A3Ptr[i];
//	}

	edge2VertexIteration_host(edgeData, vertexData, operatorFunPointer_host,
			   	operatorFunPointer_slave, &MLBParas);

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
		(this->getTopology()->getVertexNumber()*sizeof(swFloat));
	for(iter = this->getVertexMap().begin();
				iter!=this->getVertexMap().end();iter++)
	{
		vertexData->A4Ptr[iter->first]=vertexData->A1Ptr[iter->second];
	}
}

void MultiLevelBlockIterator::vertex2EdgeIteration(Arrays* neighbourData,
			Arrays* vertexData,void(*operatorFunPointer_host)
			(MLBFunParameters *MLBFunParas), void(*operatorFunPointer_slave)
			(MLBFunParameters *MLBFunParas))
{
	MLBParameters MLBParas;
	MLBParas.blockStarts       = this->getBlockStartsV();
	MLBParas.vertexStarts      = this->getVertexStarts();
	MLBParas.firstEdgeVertices
		= this->getTopology()->getFirstEdgeVertices();
	MLBParas.vertexNeighbor
		= this->getTopology()->getVertexNeighbours();
	MLBParas.cpeBlockNum  = this->getCpeBlockNum();
	MLBParas.mshBlockNum  = this->getMshBlockNum();
	MLBParas.mtxBlockNum  = this->getMtxBlockNum();
	MLBParas.maxXNum      = this->getMaxXNum();
	MLBParas.maxCells     = this->getMaxCells();
	MLBParas.maxEdges     = this->getMaxEdges();

	vertex2EdgeIteration_host(neighbourData, vertexData,
				operatorFunPointer_host, operatorFunPointer_slave,
				&MLBParas);

	map<swInt, swInt>::iterator iter;
	vertexData->A4Ptr=(swFloat*)malloc
		(this->getTopology()->getVertexNumber()*sizeof(swFloat));
	for(iter = this->getVertexMap().begin();
				iter!=this->getVertexMap().end();iter++)
	{
		if(iter->first==0) printf("%d\n",iter->second);
		vertexData->A4Ptr[iter->first]=vertexData->A1Ptr[iter->second];
	}
	printf("%f\n",vertexData->A4Ptr[0]);

	LOG("!!!---Not implemented in MLB");
//	operatorFunPointer(edgeData, vertexData, this->_topo);
}

void MultiLevelBlockIterator::reorderEdgeDataUnsymm(Arrays* edgeData)
{

}

void MultiLevelBlockIterator::reorderEdgeData(Arrays* edgeData)
{
	swFloat* lower = edgeData->A1Ptr;
	swFloat* upper = edgeData->A2Ptr;
	LOG("reorderEdgeData");
	map<swInt, swInt>::iterator iter;
	swFloat* tmpL=(swFloat*)malloc
		(this->getTopology()->getEdgeNumber()*sizeof(swFloat));
	swFloat* tmpU=(swFloat*)malloc
		(this->getTopology()->getEdgeNumber()*sizeof(swFloat));
	for(iter = this->getEdgeMap().begin();
				iter!=this->getEdgeMap().end();iter++)
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
	for(int i=0;i<this->getTopology()->getEdgeNumber();i++)
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
	if(vertexData->A1Ptr!=NULL) reorderVertexArray(vertexData->A1Ptr);
	if(vertexData->A2Ptr!=NULL) reorderVertexArray(vertexData->A2Ptr);
	if(vertexData->A3Ptr!=NULL) reorderVertexArray(vertexData->A3Ptr);
	if(vertexData->A4Ptr!=NULL) reorderVertexArray(vertexData->A4Ptr);
}


/*#########################################################
################## HERE IS THE HEADER ######################
Direct segment describes a decomposing method that decompose
the unstructure data set directly with the original
sequence the data in. 
Two level irregular tile was used implicitly to meet load
balance and efficient cache usage.
Data swapping between parallel threads was used to simulate
shared cache.
Minimize the topology matrix bandwidth or block the segments 
with high cohesion to ensure high efficiency of the Iterator 
##########################################################*/

#include <vector>
#include <sys/time.h>
#include <assert.h>
#include <climits>
#include "swMacro.h"
#include "iterator.h"
#include "directSegmentIterator.hpp"
#include "directSegmentIterator.h"
#include "register.H"

extern "C"
{
#include <athread.h>
void slave_directSegmentIterator_e2v_slave(DS_edge2VertexPara *para);
}

namespace UNAT
{

//******************************** Constructor ********************************/
DirectSegmentIterator::DirectSegmentIterator(Topology &topo, 
			swInt* vertexWeights,
			swInt* edgeWeights)
	: Iterator(topo, vertexWeights, edgeWeights)
{
	// gather statistics data
	Topology* topo_ = getTopology();
	swInt vertexNumber = topo_->getVertexNumber();
	swInt edgeNumber = topo_->getEdgeNumber();
	swInt* startVertex = topo_->getStartVertices();
	swInt* endVertex = topo_->getEndVertices();
printf("edgeNumber is %d\n", edgeNumber);
printf("vertexNumber is %d\n", vertexNumber);
//printArray("%d", startVertex, edgeNumber);
//printArray("%d", endVertex, edgeNumber);
//printArray("%d", vertexWeights, vertexNumber);
//printArray("%d", edgeWeights, edgeNumber);

	// add up the data weights on vertex and conneted edges	
	swInt dataWeights[vertexNumber];
	memset(dataWeights, 0, vertexNumber*sizeof(swInt));
	for(swInt iedge = 0; iedge<edgeNumber; iedge++)
	{
		dataWeights[startVertex[iedge]] += edgeWeights[iedge]
			+ vertexWeights[endVertex[iedge]];
		dataWeights[endVertex[iedge]] += edgeWeights[iedge]
			+ vertexWeights[startVertex[iedge]];
	}

	// segmentation
	subSegNum_ = BLOCKNUM64K;
	swInt subSegDataLimit = SAFELDM/sizeof(swInt);
//printArray("%d", dataWeights, vertexNumber);
	decomposeArray(subSegNum_, subSegDataLimit, dataWeights, vertexNumber);
	
//printArray("%d", segStarts_, segNum_*subSegNum_+1);
	// decompose edge
	edgeStarts_ = new swInt[segNum_*subSegNum_+1];
	edgeStarts_[0] = 0;
	for(swInt istart = 1; istart < segNum_*subSegNum_+1; istart++)
	{
		swInt posi;
		if(segStarts_[istart] < vertexNumber)
		{	  
//printf("*searching %d, from %dth edge, a %d length array\n",
//		   	segStarts_[istart]-1, edgeStarts_[istart-1], edgeNumber-edgeStarts_[istart-1]);
			
			// take care of the edge starts jump over segment starts
			if(startVertex[edgeStarts_[istart-1]] > segStarts_[istart]-1)
			{
#ifdef DEBUG
				printf("\n***Warning: edge starts jump over segment starts!\n\n");
#endif
				edgeStarts_[istart] = segStarts_[istart-1];
				break;
			}

			biSearch(posi, &startVertex[edgeStarts_[istart-1]], 
						segStarts_[istart]-1, edgeNumber-edgeStarts_[istart-1]);
			
			swInt findValue = segStarts_[istart]-1;
			while(posi == -1)
			{
#ifdef DEBUG
				printf("\n***Warning: vertex %d is not start vertex, findValue decreased by 1!\n\n",
							findValue);
#endif
				findValue--;
				biSearch(posi, &startVertex[edgeStarts_[istart-1]], 
							findValue, edgeNumber-edgeStarts_[istart-1]);	
			}
//printf("*at %d\n", posi);
			posi += edgeStarts_[istart-1];
			if( posi < edgeStarts_[istart-1] )
			{
				dumpError("Wrong Edge starts found!");
				exit(-1);
			}
			while( posi < edgeNumber-1
						&& startVertex[posi] == startVertex[posi+1])
				posi++;
			edgeStarts_[istart] = posi+1;
		}
		else
			edgeStarts_[istart] = edgeNumber;
	}
//printArray("%d", edgeStarts_, segNum_*subSegNum_+1);

	// construct the segment connectivity topology
	edgeNeiSeg_ = new swInt[edgeNumber];
	std::vector<swInt> segEdgeNumTmp;
	std::vector<swInt> startSeg, endSeg;
	swInt connectCount = 0;
	for(swInt iseg = 0; iseg < segNum_*subSegNum_; iseg++)
	{
		// TODO
		// Assume the count of sparse segment won't exceed 64;
//		swInt *sparSegMap = new swInt(subSegNum_);
//		swInt sparSegIdx;
//		sparSegIdx = subSegNum_;
//		memset(sparSegMap,subSegNum_,subSegNum_);
		bool hasSparseEdge = false;
		// loop sub-segment to calculate edge number between sub-segments
		std::vector<swInt> localSegEdgeNumTmp(segNum_*subSegNum_, 0);
//printf("in segrow %d(%d Edge): ", iseg, edgeStarts_[iseg+1]-edgeStarts_[iseg]);
		for(swInt iedge = edgeStarts_[iseg]; iedge < edgeStarts_[iseg+1]; iedge++)
		{
			swInt endV = endVertex[iedge];

			swInt ubound, lbound, mbound;
			lbound = 0;
			ubound = segNum_*subSegNum_;
			while(1)
			{
				if(lbound == ubound - 1) 
					break;
				mbound = (lbound+ubound)/2;
				if(segStarts_[mbound] > endV)
				  ubound = mbound;
				else if(segStarts_[mbound] <= endV)
				  lbound = mbound;
			}
			localSegEdgeNumTmp[lbound]++;
//			assert(lbound>=iseg*subSegNum_);
			if(lbound < (iseg/subSegNum_+1)*subSegNum_)
			{
				edgeNeiSeg_[iedge] = lbound - (iseg/subSegNum_)*subSegNum_;
			} else
			{
				edgeNeiSeg_[iedge] = subSegNum_;
			}
//printf("(%d,%d); ", startVertex[iedge], endV);
		}
//printf("\n\t segEdges: ");
		// compress sub-segments connection info into coo
		for(swInt jseg=0; jseg < segNum_*subSegNum_; jseg++)
		{
//if(localSegEdgeNumTmp[jseg] != 0 && jseg == iseg)
//printf("-(%d, %d); ", jseg, localSegEdgeNumTmp[jseg]);
			if(localSegEdgeNumTmp[jseg] != 0
	  					&& jseg != iseg)
			{
//printf("(%d, %d); ", jseg, localSegEdgeNumTmp[jseg]);
				segEdgeNumTmp.push_back(localSegEdgeNumTmp[jseg]);
				startSeg.push_back(iseg);
				endSeg.push_back(jseg);
				connectCount++;
			}
		}
//printf("\n");
	}
//	printArray("%d",edgeNeiSeg_,edgeNumber);

	// store the connectivity	
	segEdgeNum_ = new swInt[connectCount];
	memcpy(segEdgeNum_, &segEdgeNumTmp[0], sizeof(swInt)*connectCount);
	segConnetion_ = Topology::constructFromEdge(&startSeg[0], &endSeg[0],
				   	connectCount);

	// get the maximum edges of sparse block  in every segment 
	int maxRowEdges=0;
	this->maxRowEdges_ = 0;
	this->wrtStarts_
		= (swInt*)malloc((this->segNum_*this->subSegNum_+1)*sizeof(swInt));
	this->recvStarts_
		= (swInt*)malloc((this->segNum_*this->subSegNum_+1)*sizeof(swInt));
	this->wrtStarts_[0] = 0;
	this->recvStarts_[0] = 0;
	int maxEdges=0;
	int minEdges=10000;
	for(int iseg=0;iseg<this->segNum_;iseg++)
	{
		for(int isubseg=0;isubseg<this->subSegNum_;isubseg++)
		{
			this->recvStarts_[iseg*BLOCKNUM64K+isubseg+1]=0;
		}
	}
	for(int iseg=0;iseg<this->segNum_;iseg++)
	{
		for(int isubseg=0;isubseg<this->subSegNum_;isubseg++)
		{
			int edgeLen = this->edgeStarts_[iseg*BLOCKNUM64K+isubseg+1]
				- this->edgeStarts_[iseg*BLOCKNUM64K+isubseg];
			maxEdges = maxEdges > edgeLen ? maxEdges : edgeLen;
			minEdges = minEdges < edgeLen ? minEdges : edgeLen;
			int blockIdx = (iseg*BLOCKNUM64K+isubseg)/this->subSegNum_;
			for(int icol=this->segConnetion_->getAccuStartVertexNumbers()[iseg*BLOCKNUM64K+isubseg];icol<this->segConnetion_->getAccuStartVertexNumbers()[iseg*BLOCKNUM64K+isubseg+1];icol++)
			{
				int colIdx = this->segConnetion_->getEndVertices()[icol];
				if(colIdx >= BLOCKNUM64K*(blockIdx+1))
				{
					maxRowEdges += this->segEdgeNum_[icol];
				} else
				{
					this->recvStarts_[colIdx] += this->segEdgeNum_[icol];
				}
			}
			this->wrtStarts_[iseg*BLOCKNUM64K+isubseg+1]
				= this->wrtStarts_[iseg*BLOCKNUM64K+isubseg]
				+ maxRowEdges;
//printf("iseg: %d, isubseg: %d, maxRowEdges: %d, wrtStarts: %d\n",iseg, isubseg, maxRowEdges, this->wrtStarts_[iseg*BLOCKNUM64K+isubseg+1]);
			maxRowEdges = 0;
		}
		int segLen = this->wrtStarts_[(iseg+1)*BLOCKNUM64K]
			- this->wrtStarts_[iseg*BLOCKNUM64K];
		this->maxRowEdges_
			= this->maxRowEdges_ > segLen ? this->maxRowEdges_ : segLen;
	}
//	for(int iseg=0;iseg<this->segNum_;iseg++)
//	{
//		for(int isubseg=0;isubseg<this->subSegNum_;isubseg++)
//		{
//			printf("%d, %d\n",iseg*BLOCKNUM64K+isubseg+1, this->recvStarts_[iseg*BLOCKNUM64K+isubseg+1]);
//		}
//	}

	// Data for computing in master core
	swInt64 sparseEdgeNum = this->wrtStarts_[this->segNum_*BLOCKNUM64K];
	this->rOwner_ = (swInt*)malloc(sparseEdgeNum*sizeof(swInt));
	this->rNeighbor_ = (swInt*)malloc(sparseEdgeNum*sizeof(swInt));

	this->accuColNum_ = (swInt*)malloc(this->subSegNum_*this->segNum_*(this->subSegNum_+1)*sizeof(swInt));
	swInt *colNum = (swInt*)malloc(this->subSegNum_*this->segNum_*this->subSegNum_*sizeof(swInt));
	int segNum = this->segNum_*this->subSegNum_;
//	printf("segNum: %d\n",segNum);
	for(int iseg=0;iseg<segNum;iseg++)
	{
		int blockIdx = iseg/BLOCKNUM64K;
		int edgeNum = this->edgeStarts_[iseg+1]-this->edgeStarts_[iseg];
		for(int icol=0;icol<BLOCKNUM64K;icol++)
		{
			int idx = iseg*this->subSegNum_+icol;
			colNum[idx] = 0;
		}
		for(int icol=this->segConnetion_->getAccuStartVertexNumbers()[iseg];icol<this->segConnetion_->getAccuStartVertexNumbers()[iseg+1];icol++)
		{
			int colIdx = this->segConnetion_->getEndVertices()[icol];
			edgeNum -= this->segEdgeNum_[icol];
			if(colIdx<(blockIdx+1)*BLOCKNUM64K && colIdx>=blockIdx*BLOCKNUM64K)
			{
				int idx = iseg*this->subSegNum_+colIdx-blockIdx*BLOCKNUM64K;
				colNum[idx] = this->segEdgeNum_[icol];
//if(iseg==64) printf("%d,%d,%d,%d,%d\n",iseg,icol,colIdx,blockIdx,colNum[idx]);
			}
		}
		int idx = iseg*this->subSegNum_+iseg-blockIdx*BLOCKNUM64K;
		colNum[idx] = edgeNum;
	}
	for(int iseg=0;iseg<this->segNum_*this->subSegNum_;iseg++)
	{
		this->accuColNum_[iseg*(this->subSegNum_+1)]=0;
		for(int icol=0;icol<this->subSegNum_;icol++)
		{ 
			int idx = iseg*(this->subSegNum_+1)+icol;
			this->accuColNum_[idx+1]
				= this->accuColNum_[idx]+colNum[idx-iseg];
//			if(iseg==53) printf("%d,%d,%d,%d\n",iseg,icol,this->accuColNum_[idx+1],colNum[idx]);
		}
	}


	// Init the send/recv topology in register communication
	initOwnNeiSendList();
	initSendList(ownNeiSendList, dataList, segNum_);
//printArray("%d", segEdgeNum_, connectCount);
}

// Copy constructors
DirectSegmentIterator::DirectSegmentIterator(const DirectSegmentIterator & iterator)
{
	// To do
	dumpError("To Be Implemented!");
	exit(-1);
}

DirectSegmentIterator& DirectSegmentIterator::operator=(
			const DirectSegmentIterator& )
{
	//To do
	dumpError("To Be Implemented!");
	exit(-1);
}

//******************************** deconstructor ********************************/
DirectSegmentIterator::~DirectSegmentIterator()
{
	delete segStarts_;
	delete edgeStarts_;
	delete segEdgeNum_;
	delete edgeNeiSeg_;
	delete segEdgeNum_;
	delete segConnetion_;
}	
//******************************* Iteration processor ************************/
// Edge to vertex iterator
void DirectSegmentIterator::edge2VertexIteration(Arrays* backEdgeData, 
			Arrays* frontEdgeData, Arrays* selfConnData, Arrays* vertexData, 
			e2v_hostFunPtr fun_host, e2v_slaveFunPtr fun_slave)
{
	// To do - compute data size and redecompose segments
	// ... ...
	
	// Data for computing in master core
	this->rBackEdgeData_  = new Arrays();
	this->rFrontEdgeData_ = new Arrays();
	constructFromArrays(backEdgeData,  this->rBackEdgeData_,
				this->maxRowEdges_);
	constructFromArrays(frontEdgeData, this->rFrontEdgeData_,
				this->maxRowEdges_);
	
	// collect data
	DS_edge2VertexPara para = 
	{
		// topology relevent
		this->getTopology()->getStartVertices(),
		this->getTopology()->getEndVertices(),
		this->segNum_,
		this->subSegNum_,
		this->segStarts_,
		this->edgeStarts_,

		this->segConnetion_->getStartVertices(),
		this->segConnetion_->getEndVertices(),
		this->segConnetion_->getAccuStartVertexNumbers(),
		this->segEdgeNum_,
		this->edgeNeiSeg_,
		this->accuColNum_,
		this->wrtStarts_,
		this->recvStarts_,
		this->getTopology()->getVertexNumber(),

		// Run-time data
		0,
	
		// computing data
		backEdgeData,
		frontEdgeData,
		selfConnData,
		vertexData,

		// sorted data
		this->rOwner_,
		this->rNeighbor_,
		this->rBackEdgeData_,
		this->rFrontEdgeData_,
	
		// computing pointers
		fun_host,
		fun_slave
	};


	swFloat* x       = accessArray(vertexData,0);
	swFloat* b       = accessArray(vertexData,0);
	swFloat* diag    = accessArray(selfConnData,0);
	swFloat* lower   = accessArray(backEdgeData,0);
	swFloat* upper   = accessArray(frontEdgeData,0);
	swInt* owner     = this->getTopology()->getStartVertices();
	swInt* neighbor  = this->getTopology()->getEndVertices();
	swInt edgeNumber = getArraySize(frontEdgeData);
	Arrays tmpBackEdgeData;
	tmpBackEdgeData.fArraySizes = this->rBackEdgeData_->fArraySizes;
	tmpBackEdgeData.fArrayNum   = this->rBackEdgeData_->fArrayNum;
	tmpBackEdgeData.fArrayDims = this->rBackEdgeData_->fArrayDims;
	tmpBackEdgeData.fArrayInOut = this->rBackEdgeData_->fArrayInOut;
	tmpBackEdgeData.floatArrays
		= (swFloat**)malloc(tmpBackEdgeData.fArrayNum*sizeof(swFloat*));

	int minCol;
	for(int spIndex=this->segNum_-1;spIndex>=0;spIndex--)
	{
//		printf("init: %d\n",spIndex);
		initTable(spIndex);
		para.spIndex = spIndex;
//printf("segEdgeNum[0]: %d\n",this->segEdgeNum_[0]);
		__real_athread_spawn((void*)slave_directSegmentIterator_e2v_slave, &para);
		if(spIndex<this->segNum_-2)
		{
			int index=0;
			int startIdx = this->wrtStarts_[(spIndex+1)*BLOCKNUM64K];
			int endIdx   = this->wrtStarts_[(spIndex+2)*BLOCKNUM64K];
			for(int iArray=0;iArray<this->rBackEdgeData_->fArrayNum;iArray++)
			{
				int dims = getArrayDims(this->rBackEdgeData_,iArray);
				tmpBackEdgeData.floatArrays[iArray]
					= &this->rBackEdgeData_->floatArrays[iArray][0]
					+ startIdx*dims;
			}

			tmpBackEdgeData.fArraySizes  = endIdx-startIdx;
			this->rFrontEdgeData_->fArraySizes = 0;
			selfConnData->fArraySizes = 0;

			fun_host(&tmpBackEdgeData,this->rFrontEdgeData_,selfConnData,
						vertexData,&this->rOwner_[startIdx],
						&this->rNeighbor_[startIdx]);

//			swFloat* rLower = accessArray(&tmpBackEdgeData, 0);
//			swFloat* rUpper = accessArray(this->rFrontEdgeData_, 0);
//			int dims = this->rBackEdgeData_->fArrayDims[0];
////			int startIdx = this->wrtStarts_[(spIndex+1)*BLOCKNUM64K];
////			int endIdx   = this->wrtStarts_[(spIndex+2)*BLOCKNUM64K];
//
////        	this->rBackEdgeData_->fArraySizes
////        		= this->wrtStarts_[(spIndex+1)*BLOCKNUM64K]
////				- this->wrtStarts_[(spIndex)*BLOCKNUM64K];
////        	this->rFrontEdgeData_->fArraySizes
////        		= this->wrtStarts_[(spIndex+1)*BLOCKNUM64K]
////				- this->wrtStarts_[(spIndex)*BLOCKNUM64K];
////			printf("spIndex: %d, edgeNum: %d\n",spIndex+1,endIdx-startIdx);
////printArray("%d",this->rOwner_,edgeNumber);
//			for(int i=0;i<endIdx-startIdx;i++)
//			{
//				for(int iDim=0;iDim<dims;iDim++)
//				{
////					b[this->rOwner_[i]*dims+iDim]
////						+= rUpper[i*dims+iDim]
////						*  x[this->rNeighbor_[i]*dims+iDim];
//					b[this->rNeighbor_[i+startIdx]*dims+iDim]
//						+= rLower[i*dims+iDim]
//						*  x[this->rOwner_[i+startIdx]*dims+iDim];
//				}
//			}
		}
		
//printf("segEdgeNum[0]: %d\n",this->segEdgeNum_[0]);
//			int iMaster = spIndex+1;
//			minCol = this->segStarts_[(iMaster+1)*BLOCKNUM64K];
//			int dims = backEdgeData->fArrayDims[0];
//			for(int iseg=iMaster*BLOCKNUM64K;iseg<(iMaster+1)*BLOCKNUM64K;
//						iseg++)
//			{
//				for(int iedge=this->edgeStarts_[iseg];
//							iedge<this->edgeStarts_[iseg+1];iedge++)
//				{
//					if(neighbor[iedge]>=minCol)
//					{
//						for(int iDim=0;iDim<dims;iDim++)
//						{
////if(this->rOwner_[index] != owner[iedge]) printf("%d,%d,%d,%d\n",index,iedge,this->rOwner_[index],owner[iedge]);
////assert(this->rNeighbor_[index] == neighbor[iedge]);
////assert(rLower[index] == lower[iedge]);
////assert(rUpper[index] == upper[iedge]);
//
////						index++;
////				if(owner[iedge]==60) printf("owner: %d,%f,%f,%f\n",iedge,b[owner[iedge]],upper[iedge],x[neighbor[iedge]]);
////						b[owner[iedge]] += upper[iedge]*x[neighbor[iedge]];
////				if(neighbor[iedge]==60) printf("neighbor:%d,%f,%f,%f\n",iedge,b[neighbor[iedge]],lower[iedge],x[owner[iedge]]);
//							b[neighbor[iedge]*dims+iDim] += lower[iedge*dims+iDim]*x[owner[iedge]*dims+iDim];
//						}
//					}
//				}
//			}
//		}
		athread_join();
		destroyTable(spIndex);
//		printf("destroy: %d\n",spIndex);
	}
//	int index = 0;
////	minCol = this->segStarts_[BLOCKNUM64K];
	int startIdx = this->wrtStarts_[0];
	int endIdx   = this->wrtStarts_[BLOCKNUM64K];
	for(int iArray=0;iArray<this->rBackEdgeData_->fArrayNum;iArray++)
	{
		int dims = getArrayDims(this->rBackEdgeData_,iArray);
		this->rBackEdgeData_->floatArrays[iArray]
			+= startIdx*dims;
	}
	this->rBackEdgeData_->fArraySizes  = endIdx-startIdx;
	this->rFrontEdgeData_->fArraySizes = 0;
	selfConnData->fArraySizes = 0;

	fun_host(this->rBackEdgeData_,this->rFrontEdgeData_,selfConnData,
				vertexData,&this->rOwner_[startIdx],
				&this->rNeighbor_[startIdx]);

//	swFloat* rLower = accessArray(this->rBackEdgeData_, 0);
//	swFloat* rUpper = accessArray(this->rFrontEdgeData_, 0);
////	this->rBackEdgeData_->fArraySizes
////		= this->wrtStarts_[BLOCKNUM64K]-this->wrtStarts_[0];
////	this->rFrontEdgeData_->fArraySizes
////		= this->wrtStarts_[BLOCKNUM64K]-this->wrtStarts_[0];
//	int dims = this->rBackEdgeData_->fArrayDims[0];
//	int startIdx = this->wrtStarts_[0];
//	int endIdx   = this->wrtStarts_[BLOCKNUM64K];
////	edgeNumber = getArraySize(this->rBackEdgeData_);
////printArray("%d",this->rOwner_,edgeNumber);
//	for(int i=startIdx;i<endIdx;i++)
//	{
//		for(int iDim=0;iDim<dims;iDim++)
//		{
////			b[this->rOwner_[i]*dims+iDim]
////				+= rUpper[i*dims+iDim]*x[this->rNeighbor_[i]*dims+iDim];
//			b[this->rNeighbor_[i]*dims+iDim]
//				+= rLower[i*dims+iDim]*x[this->rOwner_[i]*dims+iDim];
//		}
//	}
//	for(int iseg=0;iseg<BLOCKNUM64K;iseg++)
//	{
//		for(int iedge=this->edgeStarts_[iseg];
//					iedge<this->edgeStarts_[iseg+1];iedge++)
//		{
//			if(neighbor[iedge]>=minCol)
//			{
//				assert(this->rOwner_[index]    == owner[iedge]);
//				assert(this->rNeighbor_[index] == neighbor[iedge]);
//				assert(rLower[index] == lower[iedge]);
//				assert(rUpper[index] == upper[iedge]);
////				printf("%d,%d, owner: (%d,%d), neighbor: (%d,%d)\n",index,iedge,this->rOwner_[index],owner[iedge],this->rNeighbor_[index],neighbor[iedge]);
//				index++;
////				if(owner[iedge]==60) printf("owner: %d,%f,%f,%f\n",iedge,b[owner[iedge]],upper[iedge],x[neighbor[iedge]]);
//				b[owner[iedge]] += upper[iedge]*x[neighbor[iedge]];
////				if(neighbor[iedge]==60) printf("neighbor:%d,%f,%f,%f\n",iedge,b[neighbor[iedge]],lower[iedge],x[owner[iedge]]);
//				b[neighbor[iedge]] += lower[iedge]*x[owner[iedge]];
//			}
//		}
//	}

	// C luncher
	//directSegmentIterator_e2v(para);
}

// Vertex to edge iterator
void DirectSegmentIterator::vertex2EdgeIteration(Arrays* neighbourData, 
			Arrays* vertexData, v2e_hostFunPtr fun_host, v2e_slaveFunPtr fun_slave)
{
	//To do
	dumpError("To Be Implemented!");
	exit(-1);
}

//****************************** Assitant Function ***************************/
// Array two level segmentation
void DirectSegmentIterator::decomposeArray(swInt subSegNum, swInt subSegDataLimit,
			swInt* dataWeights, swInt size)
{
	swInt64 weightSum = 0;
	for(swInt ivert = 0; ivert<size; ivert++)
		weightSum += (swInt64) dataWeights[ivert];
//printArray("%d", dataWeights, size);
//printf("subSegNum is %d, subSegDataLimit is %d, weightSum is %d\n", subSegNum, subSegDataLimit, weightSum);
	// calculate initial segment count with int roundup
	swInt segDataLimite = subSegNum*subSegDataLimit;
#ifdef DEBUG
	printf("subSegDataLimit is %d\n", subSegDataLimit);
	printf("weightSum is %d\n", weightSum);
#endif
	this->segNum_ = (weightSum + segDataLimite - 1)
		/ segDataLimite;
	if(segNum_ <= 0)
	{
		dumpError("too few vetex or too small %s",
					" total weights to do Array decomposing");
		exit(-1);
	}

	while(true)
	{
		// calculate segment starts
		swInt subSegCount = 0;
		swInt vertexCount = 0;
		swInt64 remainWeight = weightSum;
		swInt remainSegNum = segNum_;
		swFloat subSegVol;
		// initialize temporary segment starts 
		swInt segStartsTmp[segNum_*subSegNum+1];
		segStartsTmp[0] = 0;
		// loop segment, calculate dynamic sub-segment volume 
		while(remainSegNum > 0)
		{
			subSegVol = (swFloat)remainWeight /
				( (swFloat)remainSegNum * (swFloat)(subSegNum) );
			if(subSegVol > subSegDataLimit)
				subSegVol = subSegDataLimit;
#ifdef DEBUG	
			printf("subSegVol is %f\n", subSegVol);
#endif
			// loop the leaf sub-segment in one segment to calculate starts
			for(swInt isub = subSegCount; isub < subSegCount+subSegNum; isub++)
			{
				subSegVol = (swFloat)remainWeight /
					( (swFloat)remainSegNum * (swFloat)(subSegNum) 
					  - (swFloat) isub + (swFloat) subSegCount );

				swInt accmltWeight = 0;
				// loop vertex to accumulate sub-segment weights
				while(accmltWeight < subSegVol && vertexCount < size)
				{
					remainWeight -= dataWeights[vertexCount];
					accmltWeight += dataWeights[vertexCount];
					vertexCount++;
				}
				// check last aggressive(greedy) step
				if( accmltWeight >= subSegDataLimit )
				{
					remainWeight += dataWeights[vertexCount];
					accmltWeight -= dataWeights[vertexCount];
					vertexCount--;
				}
				segStartsTmp[isub+1] = vertexCount;
			}	
			subSegCount += subSegNum;
			remainSegNum--;
		}
//printArray("%d", segStartsTmp, segNum_*subSegNum+1);
		// check the completeness of the decomposing
#ifdef DEBUG
		printf("segNum_ is %d\n", segNum_);
#endif
		if(vertexCount < size)
		{
			swInt lastSubWeight = 0;
			swInt lastSegStart = segStartsTmp[segNum_*subSegNum-1];
			for(swInt ivert = lastSegStart; ivert < size; ivert++)
				lastSubWeight += dataWeights[ivert];
			if(lastSubWeight < subSegDataLimit)
			{
				segStartsTmp[segNum_*subSegNum] = size;
				swInt startsNum = segNum_*subSegNum+1;
				segStarts_ = new swInt[startsNum];
				memcpy(segStarts_, &segStartsTmp[0], sizeof(swInt)*startsNum);
				break;
			}
			else
				segNum_++;
		}
		else
		{
			swInt startsNum = segNum_*subSegNum+1;
			segStarts_ = new swInt[startsNum];
			memcpy(segStarts_, &segStartsTmp[0], sizeof(swInt)*startsNum);
			break;
		}
	}
}

// Init the send/receive topology during the register communication
void DirectSegmentIterator::initOwnNeiSendList()
{
	swInt segNum;
	segNum    = this->segNum_;

	swInt *owner,*neighbor,*segStarts,*edgeStarts;
	owner      = this->getTopology()->getStartVertices();
	neighbor   = this->getTopology()->getEndVertices();
	segStarts  = this->segStarts_;
	edgeStarts = this->edgeStarts_;
	
	dataList.resize(segNum*BLOCKNUM64K);
	ownNeiSendList
		= (int*)malloc(segNum*BLOCKNUM64K*BLOCKNUM64K*sizeof(int));
	for(int i=0;i<segNum*BLOCKNUM64K*BLOCKNUM64K;i++){ownNeiSendList[i]=0;}

	// Get the count of edges of every column block.
	for(int iseg=0;iseg<segNum*BLOCKNUM64K;iseg++)
	{
		int blockIdx = iseg/BLOCKNUM64K;
//		printf("%d,%d,%d\n",iseg,edgeStarts[iseg],edgeStarts[iseg+1]);
		for(int iedge=edgeStarts[iseg];iedge<edgeStarts[iseg+1];iedge++)
		{
			for(int col=blockIdx*BLOCKNUM64K;
						col<(blockIdx+1)*BLOCKNUM64K;col++)
			{
				if(neighbor[iedge]<segStarts[col+1]
							&& neighbor[iedge]>=segStarts[col])
				{
					int idx = iseg*BLOCKNUM64K+col-blockIdx*BLOCKNUM64K;
					ownNeiSendList[idx]++;
				}
			}
		}
	}

	// Get the count of packages during the register communication.
	// Get the count of efficient data in every package.
	for(int iseg=0;iseg<segNum*BLOCKNUM64K;iseg++)
	{
		int blockIdx = iseg/BLOCKNUM64K;
		int irow = iseg-blockIdx*BLOCKNUM64K;
		for(int icol=irow;icol<BLOCKNUM64K;icol++)
		{
			int idx = iseg*BLOCKNUM64K+icol;
			int symmIdx = (icol+blockIdx*BLOCKNUM64K)*BLOCKNUM64K+irow;
			ownNeiSendList[symmIdx] = ownNeiSendList[idx];
		}
	}

	for(int iseg=0;iseg<segNum*BLOCKNUM64K;iseg++)
	{
		int blockIdx = iseg/BLOCKNUM64K;
		int irow = iseg-blockIdx*BLOCKNUM64K;
		for(int icol=0;icol<BLOCKNUM64K;icol++)
		{
			int idx = iseg*BLOCKNUM64K+icol;
			if(icol==irow)
			{
				ownNeiSendList[idx] = 0;
				continue;
			}
			if(ownNeiSendList[idx]==0) continue;
			int edgeNum = ownNeiSendList[idx];
			while(edgeNum>6)
			{
				dataList[iseg].push_back(6);
				edgeNum -=6;
			}
			dataList[iseg].push_back(edgeNum);
			int packNum = (ownNeiSendList[idx]-1)/6+1;
//			int symmIdx = (icol+blockIdx*BLOCKNUM64K)*BLOCKNUM64K+irow;
//			printf("%d,%d,%d,%d\n",iseg,icol,ownNeiSendList[idx],packNum);
			ownNeiSendList[idx] = packNum;
//			ownNeiSendList[symmIdx] = packNum;
		}
	}

//	for(int iseg=0;iseg<segNum*BLOCKNUM64K;iseg++)
//	{
//		int blockIdx = iseg/BLOCKNUM64K;
//		int irow = iseg-blockIdx*BLOCKNUM64K;
//		for(int icol=0;icol<BLOCKNUM64K;icol++)
//		{
//			int idx = iseg*BLOCKNUM64K+icol;
//		}
//	}
}

} //namespace UNAT


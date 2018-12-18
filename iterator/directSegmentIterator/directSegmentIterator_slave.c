#include "swMacro.h"
#include "slave.h"
#include <assert.h>
#include "rlmpi.h"
#include "slaveUtils.h"
#include "directSegmentIterator.h"
#include "iterator.h"

//// Topology
//static __thread_local_fix swInt   *segStarts,*edgeStarts,*owner,*neighbor;
//// Field data
//static __thread_local_fix Arrays *upper,*lower,*diag,*vertex;
//// spIndex: index of segment
//// maxEdges: maximum count of edges among the subsegments in same segment
//// isXExist: whether x/b is contained in vertex data.
//// dimension: the dimension of data; 1--x; 2--x,y; 3--x,y,z;
//static __thread_local_fix swInt   spIndex,maxEdges,cpeBlockNum,isXExist,dimension;
////
//static __thread_local_fix swInt   vertexNumber,edgeNumber;
//
//// variables(***_slave) are allocated in LDM
//// sVertex_slave: local vertex data;
//// rVertex_slave: global vertex data including received data from other
//// blocks for the arrays which fArratInOut[i]=COPYIN, while others keep
//// local.
static __thread_local Arrays  sVertex_slave,diag_slave;
static __thread_local Arrays  lowVertex_slave,upVertex_slave;
static __thread_local swInt   *owner_slave,*neighbor_slave,*sNeighbor_slave;
static __thread_local swInt   *accuColNum;
//// colIdx: the column block index of each edge
//// colNum: the counts of each column block
//// accuColNum: another representation of colNum
static __thread_local swInt  *accuColNum_slave;
static __thread_local swInt  *edgeNeiSeg_slave;
// sLower_slave: lower data stored upper triangle.
// rLower_slave: lower data received from other blocks
static __thread_local Arrays  sLower_slave,upper_slave,rLower_slave;
//// record the base in quicksort
static __thread_local swInt   pivOwner,pivNeighbor;
static __thread_local swFloat **pivLower,**pivUpper;
static __thread_local volatile int myId,recvNum;

//inline void copyArray(Arrays *slave, Arrays *host, swInt length,
//			swInt startIdx, char* ldm_space, char* ldm_space_end,
//			int ldm_size);
inline void quickSort(int *arr, int elements);
inline int record(const int L);
inline int exchange(const int L, const int R);
inline int update(const int L);
//inline void allocRecvArray(Arrays *send, Arrays *recv, 
//			swInt diagNum, swInt offDiagNum, char* ldm_space,
//			char* ldm_space_ned, int ldm_size);
inline void outputArray(Arrays *slave, Arrays *host,
			swInt startIdx, swInt length);

#define copyArray(slave, host, length, startIdx) \
{ \
	int i; \
	swInt fArrayNum    = host->fArrayNum; \
	slave.fArrayNum   = fArrayNum; \
	slave.fArraySizes = length; \
	LDM_NEW(slave.fArrayDims, swInt, fArrayNum); \
	LDM_NEW(slave.fArrayInOut, swInt, fArrayNum); \
	DMA_Get(slave.fArrayDims,host->fArrayDims,fArrayNum*sizeof(swInt)); \
	DMA_Get(slave.fArrayInOut,host->fArrayInOut,fArrayNum*sizeof(swInt));\
	LDM_NEW(slave.floatArrays, swFloat*, fArrayNum); \
	for(i=0;i<fArrayNum;i++) \
	{ \
		LDM_NEW(slave.floatArrays[i], swFloat, \
					length*slave.fArrayDims[i]); \
		DMA_Get(slave.floatArrays[i], \
			&host->floatArrays[i][slave.fArrayDims[i]*startIdx], \
			length*slave.fArrayDims[i]*sizeof(swFloat)); \
	} \
}

#define allocRecvArray(recv, send, diagNum, offDiagNum)\
{ \
	int i,j; \
	recv.fArrayNum   = send.fArrayNum; \
	recv.fArraySizes = diagNum+offDiagNum; \
	LDM_NEW(recv.fArrayDims,swInt,recv.fArrayNum); \
	LDM_NEW(recv.fArrayInOut,swInt,recv.fArrayNum);\
	LDM_NEW(recv.floatArrays,swFloat*,recv.fArrayNum); \
	for(i=0;i<recv.fArrayNum;i++)\
	{\
		recv.fArrayDims[i] = send.fArrayDims[i];\
		recv.fArrayInOut[i] = send.fArrayInOut[i];\
		LDM_NEW(recv.floatArrays[i],swFloat,\
					(diagNum+offDiagNum)*recv.fArrayDims[i]);\
		for(j=0;j<diagNum*recv.fArrayDims[i];j++)\
		{\
			recv.floatArrays[i][j]\
				= send.floatArrays[i][j];\
		}\
	}\
}
void directSegmentIterator_e2v_slave(DS_edge2VertexPara *para)
{
// **********************************************************************
// Start Initiation
// TODO
// move to another function\
// **********************************************************************
	myId = athread_get_id(-1);

	INIT_LDM_SPACE(56000);
	swInt *segStarts,*edgeStarts,*owner,*neighbor,*edgeNeiSeg,*accuColNum;
	swInt *rOwner,*rNeighbor;
	segStarts  = para->segStarts;
	edgeStarts = para->edgeStarts;
	edgeNeiSeg = para->edgeNeiSeg;
	accuColNum = para->accuColNum;
	owner      = para->startVertex;
	neighbor   = para->endVertex;
	rOwner     = para->rOwner;
	rNeighbor  = para->rNeighbor;

	Arrays *upper,*lower,*diag,*vertex,*rLower,*rUpper;
	upper  = para->frontEdgeData;
	lower  = para->backEdgeData;
	diag   = para->selfConnData;
	vertex = para->vertexData;
	rLower = para->rBackEdgeData;
	rUpper = para->rFrontEdgeData;

	swInt spIndex,vertexNumber;
	spIndex = para->spIndex;
	vertexNumber = para->vertexNumber;

	volatile swInt cellLen,edgeLen;
	int i,j;
	swInt segStarts_slave[BLOCKNUM64K+1],edgeStarts_slave[BLOCKNUM64K+1];
	volatile swInt segStart = segStarts_slave[myId];
	volatile swInt segEnd   = segStarts_slave[myId+1];
	DMA_Get(&segStarts_slave[0],&segStarts[spIndex*BLOCKNUM64K],
				(BLOCKNUM64K+1)*sizeof(swInt));
	DMA_Get(&edgeStarts_slave[0],&edgeStarts[spIndex*BLOCKNUM64K],
				(BLOCKNUM64K+1)*sizeof(swInt));
//	LDM_NEW(colNum,swInt,BLOCKNUM64K+1);
	LDM_NEW(accuColNum_slave,swInt,BLOCKNUM64K+2);
	cellLen    = segStarts_slave[myId+1]-segStarts_slave[myId];

    // Get the maximum edge count of all segments.
	// assert(maxEdges>recvNum+sendNum);
    swInt maxEdges = 0;
	for(i=0;i<BLOCKNUM64K;i++)
	{
		edgeLen  = edgeStarts_slave[i+1]-edgeStarts_slave[i];
		maxEdges = maxEdges > edgeLen ? maxEdges : edgeLen;
	}
	maxEdges = maxEdges*1.5;

    edgeLen        = edgeStarts_slave[myId+1]-edgeStarts_slave[myId];
	LDM_NEW(owner_slave,swInt,maxEdges);
	LDM_NEW(neighbor_slave,swInt,maxEdges);
	LDM_NEW(sNeighbor_slave,swInt,maxEdges);
	LDM_NEW(edgeNeiSeg_slave,swInt,maxEdges);
    DMA_Get(sNeighbor_slave,&neighbor[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
    DMA_Get(owner_slave,&owner[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
    DMA_Get(edgeNeiSeg_slave,&edgeNeiSeg[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));

//    for(i=0;i<BLOCKNUM64K+1;i++) {colNum[i]=0;}
	for(i=0;i<edgeLen;i++)
	{
		if(sNeighbor_slave[i]<segStarts_slave[BLOCKNUM64K])
		{
			edgeNeiSeg_slave[i] -=spIndex*BLOCKNUM64K;
//			for(j=myId;j<BLOCKNUM64K;j++)
//			{
//				if(sNeighbor_slave[i]<segStarts_slave[j+1]
//							&& sNeighbor_slave[i]>=segStarts_slave[j])
//				{
//					colIdx[i] = j;
//					assert(edgeNeiSeg_slave[i]==j);
//					colNum[j]++;
//				}
//			}
		} else
		{
			edgeNeiSeg_slave[i]=BLOCKNUM64K;
//			colIdx[i]=BLOCKNUM64K;
//			colNum[BLOCKNUM64K]++;
		}
	}

//    for(i=0;i<BLOCKNUM64K+2;i++) {accuColNum_slave[0]=0;}
//	for(i=1;i<BLOCKNUM64K+1;i++)
//	{
//		accuColNum_slave[i+1] += accuColNum_slave[i]+colNum[i-1];
//	}
#ifdef DEBUG
//	assert(accuColNum[myId+1]==colNum[myId]);
//	assert(accuColNum[BLOCKNUM64K+1]==edgeLen);
#endif

    // Get the global index of each edge for quicksort
	// priority: column block index > row index > col index
//    for(i=0;i<edgeLen;i++)
//	{
//		edgeGlobalIdx[i]
//			= colIdx[i]*cellLen*vertexNumber
//			+ owner_slave[i]*vertexNumber
//			+ sNeighbor_slave[i];
//}
// ***********************************************************************
// END Initiation
// ***********************************************************************

	// data stored in slave core
//	Arrays sVertex_slave,lowVertex_slave,upVertex_slave,upper_slave;
//	Arrays sLower_slave,rLower_slave,diag_slave;
//	if(myId==0) printf("vertex-fArrayNum: %d\n",vertex.fArrayNum);
    copyArray(sVertex_slave, para->vertexData,
				cellLen, segStarts_slave[myId]);
	copyArray(upper_slave,   para->frontEdgeData,
				edgeLen, edgeStarts_slave[myId]);
	copyArray(sLower_slave,  para->backEdgeData,
				edgeLen, edgeStarts_slave[myId]);
	copyArray(diag_slave,    para->selfConnData,
				cellLen, segStarts_slave[myId]);

//	swInt maxDim = 0;
//	for(i=0;i<sLower_slave.fArrayNum;i++)
//	{
//		maxDim = maxDim > sLower_slave.fArrayDims[i]
//			? maxDim : sLower_slave.fArrayDims[i];
//	}
//	swInt pivLower[sLower_slave.fArrayNum][maxDim];
//	maxDim = 0;
//	for(i=0;i<upper_slave.fArrayNum;i++)
//	{
//		maxDim = maxDim > upper_slave.fArrayDims[i]
//			? maxDim : upper_slave.fArrayDims[i];
//	}
//	swInt pivUpper[upper_slave.fArrayNum][maxDim];

//	// Record the compared base in quicksort
//	LDM_NEW(pivLower,swFloat*,sLower_slave.fArrayNum);
//	LDM_NEW(pivUpper,swFloat*,upper_slave.fArrayNum);
//	for(i=0;i<sLower_slave.fArrayNum;i++)
//	{
//		LDM_NEW(pivLower[i],swFloat,sLower_slave.fArrayDims[i]);
//	}
//	for(i=0;i<upper_slave.fArrayNum;i++)
//	{
//		LDM_NEW(pivUpper[i],swFloat,upper_slave.fArrayDims[i]);
//	}

//if(myId==0) printArray("%d",sNeighbor_slave,edgeLen);

	swInt *tmpInt,*bonus;
	swFloat *tmpFloat;
	LDM_NEW(tmpInt,swInt,edgeLen);
	DMA_Get(&accuColNum_slave[1],&accuColNum[spIndex*BLOCKNUM64K*(BLOCKNUM64K+1)+myId*(BLOCKNUM64K+1)],(BLOCKNUM64K+1)*sizeof(swInt));
//	if(myId==62)
//	{
//		for(i=0;i<BLOCKNUM64K+2;i++)
//		  printf("%d\n",accuColNum_slave[i]);
//	}

	LDM_NEW(bonus,int,BLOCKNUM64K+1);

	int k,dims,idim;
	// TODO
	// Get the maximum dimension
	dims = sLower_slave.fArrayDims[0];
	LDM_NEW(tmpFloat,swFloat,edgeLen*dims);
	// Reorder lower
	for(i=0;i<sLower_slave.fArrayNum;i++)
	{
		int dims;
		dims = sLower_slave.fArrayDims[i];
		for(j=0;j<BLOCKNUM64K+1;j++) {bonus[j]=0;}
		DMA_Get(tmpFloat,&para->backEdgeData->floatArrays[i][dims*edgeStarts_slave[myId]],edgeLen*dims*sizeof(swFloat));
		for(j=0;j<edgeLen;j++)
		{
			for(idim=0;idim<dims;idim++)
			{
				sLower_slave.floatArrays[i][accuColNum_slave[edgeNeiSeg_slave[j]+1]+bonus[edgeNeiSeg_slave[j]]] = tmpFloat[j*dims+idim];
				bonus[edgeNeiSeg_slave[j]]++;
			}
		}
	}
	// Reorder Upper
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		int dims;
		dims = upper_slave.fArrayDims[i];
		for(j=0;j<BLOCKNUM64K+1;j++) {bonus[j]=0;}
		DMA_Get(tmpFloat,&para->frontEdgeData->floatArrays[i][dims*edgeStarts_slave[myId]],edgeLen*dims*sizeof(swFloat));
		for(j=0;j<edgeLen;j++)
		{
			for(idim=0;idim<dims;idim++)
			{
				upper_slave.floatArrays[i][accuColNum_slave[edgeNeiSeg_slave[j]+1]+bonus[edgeNeiSeg_slave[j]]] = tmpFloat[j*dims+idim];
				bonus[edgeNeiSeg_slave[j]]++;
			}
		}
	}

	// Reorder neighbor
	for(i=0;i<BLOCKNUM64K+1;i++) {bonus[i]=0;}
    DMA_Get(tmpInt,&neighbor[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
	for(i=0;i<edgeLen;i++)
	{
		sNeighbor_slave[accuColNum_slave[edgeNeiSeg_slave[i]+1]+bonus[edgeNeiSeg_slave[i]]] = tmpInt[i];
		bonus[edgeNeiSeg_slave[i]]++;
	}

	// Reorder owner
	for(i=0;i<BLOCKNUM64K+1;i++) {bonus[i]=0;}
    DMA_Get(tmpInt,&owner[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
	for(i=0;i<edgeLen;i++)
	{
		owner_slave[accuColNum_slave[edgeNeiSeg_slave[i]+1]] =tmpInt[i];
		accuColNum_slave[edgeNeiSeg_slave[i]+1]++;
	}


//	quickSort(&edgeGlobalIdx[0],edgeLen);

//	for(i=0;i<edgeLen;i++) {assert(test[i]==sNeighbor_slave[i]);}

	// Communicate neighbor
    swInt ownNeiSendIdx[BLOCKNUM64K],ownNeiRecvIdx[BLOCKNUM64K];
	swInt ownNeiRecvList[BLOCKNUM64K];
	int lastId,startIdx;
    for(i=0;i<accuColNum_slave[myId+1];i++)
	{
		neighbor_slave[i] = sNeighbor_slave[i];
	}
	for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
	for(i=0;i<_total_send_pcg;i++)
	{
		if(myId>_sPacks[i].dst_id) continue;
		startIdx
			=accuColNum_slave[_sPacks[i].dst_id]+ownNeiSendIdx[_sPacks[i].dst_id];
		for(j=0;j<_sPacks[i].cva;j++)
		{
			_sPacks[i].data[j] = sNeighbor_slave[startIdx+j];
		}
		ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].cva;
	}
	transform_data();
	recvNum=accuColNum_slave[myId+1];
	lastId=0;
	k=0;
	for(i=0;i<BLOCKNUM64K;i++) {ownNeiRecvList[i]=accuColNum_slave[myId+1];}
	for(i=0;i<_total_send_pcg;i++)
	{
		if(_rPacks[i].src_id>myId) continue;
		for(j=0;j<_rPacks[i].cva;j++)
		{
			neighbor_slave[recvNum] = _rPacks[i].data[j];
			recvNum++;
		}
		if(_rPacks[i].src_id!=lastId) {k++;lastId=_rPacks[i].src_id;}
		ownNeiRecvList[k+1]=recvNum;
		ownNeiRecvIdx[_rPacks[i].src_id]=k;
	}
	recvNum-=accuColNum_slave[myId+1];
//if(myId==5) printArray("%d",neighbor_slave,colNum[myId]+recvNum);
//if(myId==0) printArray("%d",colIdx,edgeLen);
#ifdef DEBUG
//	assert(maxEdges>recvNum+colNum[myId]);
#endif

	allocRecvArray(rLower_slave, sLower_slave,
				accuColNum_slave[myId+1], recvNum);
	allocRecvArray(lowVertex_slave, sVertex_slave, cellLen, recvNum);
	allocRecvArray(upVertex_slave, sVertex_slave,
				cellLen, edgeLen-accuColNum_slave[myId+1]);
	// The vertexData was seperated into two parts
	for(i=0;i<lowVertex_slave.fArrayNum;i++)
	{
		dims = lowVertex_slave.fArrayDims[i];
		if(lowVertex_slave.fArrayInOut[i]==COPYOUT)
		{
			for(j=0;j<cellLen*dims;j++)
			{
				lowVertex_slave.floatArrays[i][j] /= 2; 
				upVertex_slave.floatArrays[i][j]  /= 2;
			}
		}
	}

	// Communicate Lower
	int iArray,iDim;
    for(iArray=0;iArray<rLower_slave.fArrayNum;iArray++)
	{
		dims = rLower_slave.fArrayDims[iArray];
		for(iDim=0;iDim<dims;iDim++)
		{
			for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
			for(i=0;i<_total_send_pcg;i++)
			{
				if(myId>_sPacks[i].dst_id) continue;
				startIdx = accuColNum_slave[_sPacks[i].dst_id] * dims
					+ ownNeiSendIdx[_sPacks[i].dst_id] + iDim;
				for(j=0;j<_sPacks[i].cva;j++)
				{
					_sPacks[i].data[j]
						=sLower_slave.floatArrays[iArray][startIdx+j*dims];
				}
				ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].cva*dims;
			}
			transform_data();
			startIdx=accuColNum_slave[myId+1]*dims+iDim;
			for(i=0;i<_total_send_pcg;i++)
			{
				if(_rPacks[i].src_id>myId) continue;
				for(j=0;j<_rPacks[i].cva;j++)
				{
					rLower_slave.floatArrays[iArray][startIdx]
						= _rPacks[i].data[j];
					startIdx+=dims;
				}
			}
		}
	}
#ifdef DEBUG
//if(myId==1) printArray("%f",rLower_slave.floatArrays[0],recvNum+colNum[myId]);
//	assert(startIdx==recvNum+colNum[myId]);
#endif

	// Communicate X
	swInt idx;
	swInt ownNeiSendIdxUToL[BLOCKNUM64K],ownNeiSendIdxLToU[BLOCKNUM64K];
    for(iArray=0;iArray<sVertex_slave.fArrayNum;iArray++)
	{
		if(sVertex_slave.fArrayInOut[iArray]==COPYOUT) continue;
		dims = sVertex_slave.fArrayDims[iArray];
		for(iDim=0;iDim<dims;iDim++)
		{
			for(i=0;i<BLOCKNUM64K;i++)
			{
				ownNeiSendIdxUToL[i]=0;
				ownNeiSendIdxLToU[i]=0;
			}
			for(i=0;i<_total_send_pcg;i++)
			{
				if(myId<_sPacks[i].dst_id)
				{
					startIdx
						= accuColNum_slave[_sPacks[i].dst_id]
						+ ownNeiSendIdxUToL[_sPacks[i].dst_id];
					for(j=0;j<_sPacks[i].cva;j++)
					{
						idx=owner_slave[startIdx+j]-segStarts_slave[myId];
						_sPacks[i].data[j]=sVertex_slave.floatArrays[iArray][idx*dims+iDim];
					}
					ownNeiSendIdxUToL[_sPacks[i].dst_id]+=_sPacks[i].cva;
				} else
				{
					startIdx
						= ownNeiRecvList[ownNeiRecvIdx[_sPacks[i].dst_id]]
						+ ownNeiSendIdxLToU[_sPacks[i].dst_id];
					for(j=0;j<_sPacks[i].cva;j++)
					{
						idx=neighbor_slave[startIdx+j]-segStarts_slave[myId];
						_sPacks[i].data[j]=sVertex_slave.floatArrays[iArray][idx*dims+iDim];
//				if(myId==7 && _sPacks[i].dst_id==1) printf("packs: %d,%d,%f\n",idx,startIdx,_sPacks[i].data[j]);
					}
					ownNeiSendIdxLToU[_sPacks[i].dst_id]+=_sPacks[i].cva;
				}
			}
			transform_data();
			int upVertexIdx=cellLen*dims+iDim;
			int lowVertexIdx=cellLen*dims+iDim;
		    for(i=0;i<_total_send_pcg;i++)
			{
				if(_rPacks[i].src_id>myId)
				{
					for(j=0;j<_rPacks[i].cva;j++)
					{
//if(myId==0 && spIndex==1) printf("%d,%d,%d,%f\n",i,j,_rPacks[i].src_id,_rPacks[i].data[j]);
						upVertex_slave.floatArrays[iArray][upVertexIdx]
							= _rPacks[i].data[j];
						upVertexIdx+=dims;
					}
				} else
				{
					for(j=0;j<_rPacks[i].cva;j++)
					{
						lowVertex_slave.floatArrays[iArray][lowVertexIdx]
							= _rPacks[i].data[j];
						lowVertexIdx+=dims;
					}
				}
			}
		}
	}
//if(myId==0) printArray("%f",upVertex_slave.floatArrays[0],cellLen+edgeLen-colNum[myId]);
//if(myId==1) printArray("%d",sNeighbor_slave,edgeLen);
//if(myId==1) printArray("%f",upVertex_slave.floatArrays[0],edgeLen);

	// Transform the sorted data to master core
	swInt wrtStart = para->wrtStarts[myId+spIndex*BLOCKNUM64K];
	swInt wrtEnd   = para->wrtStarts[myId+spIndex*BLOCKNUM64K+1];
	int length = wrtEnd-wrtStart;
	// Output owner/neighbor
	if(length>0)
	{
//		if(myId==7) printf("***%d,%d,%d***\n",segStarts_slave[myId],accuColNum_slave[BLOCKNUM64K],owner_slave[accuColNum_slave[BLOCKNUM64K]]);
		DMA_Put(&rOwner[wrtStart],
				&owner_slave[accuColNum_slave[BLOCKNUM64K]],
				length*sizeof(swInt));
		DMA_Put(&rNeighbor[wrtStart],
				&sNeighbor_slave[accuColNum_slave[BLOCKNUM64K]],
				length*sizeof(swInt));
	}

	for(i=0;i<maxEdges;i++)
	{
		owner_slave[i]     -= segStarts_slave[myId];
		sNeighbor_slave[i] -= segStarts_slave[myId];
		neighbor_slave[i]  -= segStarts_slave[myId];
	}
	for(i=accuColNum_slave[myId+1];i<maxEdges;i++)
	{
		sNeighbor_slave[i] = cellLen+i-accuColNum_slave[myId+1];
	}
//if(myId==0 && spIndex==1) printArray("%d",sNeighbor_slave,edgeLen);
//if(myId==0 && spIndex==1) printArray("%f",upVertex_slave.floatArrays[0],cellLen+edgeLen-colNum[myId]);
//if(myId==1) printf("cellLen: %d, colNum[myId]: %d\n",cellLen,colNum[myId]);

	// Compute the upper triangle
    rLower_slave.fArraySizes    = 0;
	diag_slave.fArraySizes  = 0;
	upper_slave.fArraySizes = accuColNum_slave[BLOCKNUM64K];
//	para->fun_slave(&upper_slave,&rLower_slave,&diag_slave,&upVertex_slave,
//				owner_slave,sNeighbor_slave);
	//selfConn computation
	swFloat* diagA	= accessArray(&diag_slave, 0);
	swFloat* x		= accessArray(&upVertex_slave, 0);
	swFloat* b		= accessArray(&upVertex_slave, 1);
	//frontEdge computation
	swFloat* upperA	= accessArray(&upper_slave, 0);
	swInt edgeNumber = getArraySize(&upper_slave);
	swInt iedge;
	for( iedge = 0; iedge < edgeNumber; iedge++)
{
//if(owner_slave[iedge]+segStarts_slave[myId]==60) printf("upper:%d,%d,%d,%f,%f,%f\n",iedge,owner_slave[iedge],sNeighbor_slave[iedge],b[owner_slave[iedge]],upperA[iedge],x[sNeighbor_slave[iedge]]);
			b[owner_slave[iedge]] += upperA[iedge]*x[sNeighbor_slave[iedge]];
}

//if(myId==4) printArray("%f",upVertex_slave.floatArrays[1],cellLen);

	// Compute the lower triangle
	for(i=accuColNum_slave[myId+1];i<maxEdges;i++)
	{
		owner_slave[i] = cellLen+i-accuColNum_slave[myId+1];
	}
	x		= accessArray(&lowVertex_slave, 0);
	b		= accessArray(&lowVertex_slave, 1);
    upper_slave.fArraySizes  = 0;
	rLower_slave.fArraySizes = recvNum+accuColNum_slave[myId+1];
//	para->fun_slave(&upper_slave,&rLower_slave,&diag_slave,
//				&lowVertex_slave,owner_slave,neighbor_slave);
	//backEdge computation
	swFloat* lowerA	= accessArray(&rLower_slave, 0);
	edgeNumber = getArraySize(&rLower_slave);
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
//if(neighbor_slave[iedge]+segStarts_slave[myId]==60) printf("lower:%d,%d,%d,%f,%f,%f\n",iedge,owner_slave[iedge],neighbor_slave[iedge],b[neighbor_slave[iedge]],lowerA[iedge],x[owner_slave[iedge]]);
			b[neighbor_slave[iedge]] += lowerA[iedge]*x[owner_slave[iedge]];
	}
//if(myId==4) printArray("%f",lowVertex_slave.floatArrays[1],cellLen);

	// Compute the diagonal
    rLower_slave.fArraySizes = 0;
    upper_slave.fArraySizes  = 0;
	diag_slave.fArraySizes   = cellLen;
//	para->fun_slave(&upper_slave,&rLower_slave,&diag_slave,
//				&upVertex_slave,owner_slave,sNeighbor_slave);
	swInt vertexNum = getArraySize(&diag_slave);
	swInt ivertex;
	for( ivertex = 0; ivertex < vertexNum; ivertex++)
	{
//if(ivertex+segStarts_slave[myId]==60) printf("diag:%f,%f,%f\n",b[ivertex],diagA[ivertex],x[ivertex]);
		b[ivertex] += diagA[ivertex]*x[ivertex];
	}
//if(myId==4) printArray("%f",lowVertex_slave.floatArrays[1],cellLen);
	
	// Combine lower vertex and upper vertex
	for(i=0;i<lowVertex_slave.fArrayNum;i++)
	{
		dims = lowVertex_slave.fArrayDims[i];
		if(lowVertex_slave.fArrayInOut[i]==COPYOUT)
		{
			for(j=0;j<cellLen*dims;j++)
			{
				lowVertex_slave.floatArrays[i][j] 
					+= upVertex_slave.floatArrays[i][j];
			}
		}
	}

//if(myId==0) printArray("%f",lowVertex_slave.floatArrays[1],cellLen);
	// TODO
	// Assume that lower data won't be output
//if(myId==7) printf("%d,%d\n",segStart,segStarts_slave[myId]);
	lowVertex_slave.fArraySizes   = cellLen;
	diag_slave.fArraySizes        = cellLen;
	upper_slave.fArraySizes       = edgeLen;
	// Output vertexData
	for(i=0;i<lowVertex_slave.fArrayNum;i++)
	{
		DMA_Put(&(vertex->floatArrays[i]
				[segStarts_slave[myId]*lowVertex_slave.fArrayDims[i]]),
				lowVertex_slave.floatArrays[i],
				cellLen*lowVertex_slave.fArrayDims[i]*sizeof(swFloat));
	}
	// Output selfConnData
	for(i=0;i<diag_slave.fArrayNum;i++)
	{
		if(diag_slave.fArrayInOut[i]==COPYIN) continue;
		DMA_Put(&(diag->floatArrays[i]
				[segStarts_slave[myId]*diag_slave.fArrayDims[i]]),
				diag_slave.floatArrays[i],
				cellLen*diag_slave.fArrayDims[i]*sizeof(swFloat));
	}
	// Transform the sorted data to master core
	// Output upperData
	if(length<=0) return;
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		DMA_Put(&rUpper->floatArrays[i][wrtStart*upper_slave.fArrayDims[i]],
				&upper_slave.floatArrays[i]
				[accuColNum_slave[BLOCKNUM64K]*upper_slave.fArrayDims[i]],
				length*upper_slave.fArrayDims[i]*sizeof(swFloat));
	}
	// Output lowerData
	for(i=0;i<sLower_slave.fArrayNum;i++)
	{
		DMA_Put(&rLower->floatArrays[i]
				[wrtStart*sLower_slave.fArrayDims[i]],
				&sLower_slave.floatArrays[i]
				[accuColNum_slave[BLOCKNUM64K]*sLower_slave.fArrayDims[i]],
				length*sLower_slave.fArrayDims[i]*sizeof(swFloat));
	}
}




inline void quickSort(int *arr, int elements)
{
#define  MAX_LEVELS  300
	int  piv, beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap;
	beg[0]=0;
	end[0]=elements;
	while (i>=0)
	{
		L=beg[i];
		R=end[i]-1;
		if (L<R)
		{
			piv=arr[L];
			record(L);
			while (L<R)
			{
				while (arr[R]>=piv && L<R) R--;
				if (L<R) {arr[L]=arr[R]; exchange(L,R); L++;}
				while (arr[L]<=piv && L<R) L++;
				if (L<R) {arr[R]=arr[L]; exchange(R,L); R--;}
			}
			arr[L]=piv;
			update(L);
			beg[i+1]=L+1;
			end[i+1]=end[i];
			end[i++]=L;
			if (end[i]-beg[i]>end[i-1]-beg[i-1])
			{
				swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
				swap=end[i]; end[i]=end[i-1]; end[i-1]=swap;
			}
		}
		else {i--;}
	}
}

inline int record(const int L)
{
	int i,j;
//	pivOwner    = owner_slave[L];
//	pivNeighbor = sNeighbor_slave[L];
//	for(i=0;i<sLower_slave.fArrayNum;i++)
//	{
//		for(j=0;j<sLower_slave.fArrayDims[i];j++)
//		  pivLower[i][j] = sLower_slave.floatArrays[i][L+j];
//	}
//	for(i=0;i<upper_slave.fArrayNum;i++)
//	{
//		for(j=0;j<upper_slave.fArrayDims[i];j++)
//		  pivUpper[i][j] = upper_slave.floatArrays[i][L+j];
//	}
}

inline int exchange(const int L, const int R)
{
	int i,j;
//	owner_slave[L]     = owner_slave[R];
//	sNeighbor_slave[L] = sNeighbor_slave[R];
//	for(i=0;i<sLower_slave.fArrayNum;i++)
//	{
//		for(j=0;j<sLower_slave.fArrayDims[i];j++)
//		  sLower_slave.floatArrays[i][L+j]
//			  = sLower_slave.floatArrays[i][R+j];
//	}
//	for(i=0;i<upper_slave.fArrayNum;i++)
//	{
//		for(j=0;j<upper_slave.fArrayDims[i];j++)
//		  upper_slave.floatArrays[i][L+j]
//			  = upper_slave.floatArrays[i][R+j];
//	}
}

inline int update(const int L)
{
	int i,j;
//	owner_slave[L]     = pivOwner;
//	sNeighbor_slave[L] = pivNeighbor;
//	for(i=0;i<sLower_slave.fArrayNum;i++)
//	{
//		for(j=0;j<sLower_slave.fArrayDims[i];j++)
//		  sLower_slave.floatArrays[i][L+j] = pivLower[i][j];
//	}
//	for(i=0;i<upper_slave.fArrayNum;i++)
//	{
//		for(j=0;j<upper_slave.fArrayDims[i];j++)
//		  upper_slave.floatArrays[i][L+j] = pivUpper[i][j];
//	}
}



inline void outputArray(Arrays *slave, Arrays *host,
			swInt startIdx, swInt length)
{
//	int i;
//	host->fArraySizes = length;
//	for(i=0;i<slave->fArrayNum;i++)
//	{
//		if(slave->fArrayInOut[i]==COPYIN) continue;
//		DMA_Put(&(host->floatArrays[i][startIdx*slave->fArrayDims[i]]),
//				&slave->floatArrays[i][],
//				length*slave->fArrayDims[i]*sizeof(swFloat));
//	}
}


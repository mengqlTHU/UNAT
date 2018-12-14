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
static __thread_local_fix Arrays  sVertex_slave,diag_slave;
static __thread_local_fix Arrays  lowVertex_slave,upVertex_slave;
static __thread_local_fix swInt   *owner_slave,*neighbor_slave,*sNeighbor_slave;
//static __thread_local_fix swInt   *segStarts_slave,*edgeStarts_slave;
//// colIdx: the column block index of each edge
//// colNum: the counts of each column block
//// accuColNum: another representation of colNum
//static __thread_local_fix swInt  *colIdx,*colNum,*accuColNum,*edgeGlobalIdx;
// sLower_slave: lower data stored upper triangle.
// rLower_slave: lower data received from other blocks
static __thread_local_fix Arrays  sLower_slave,upper_slave,rLower_slave;
//// record the base in quicksort
static __thread_local_fix swInt   pivOwner,pivNeighbor;
static __thread_local_fix swFloat **pivLower,**pivUpper;
static __thread_local_fix volatile int myId,recvNum;

inline void copyArray(Arrays *slave, Arrays *host,
			swInt length, swInt startIdx);
inline void quickSort(int *arr, int elements);
inline int record(const int L);
inline int exchange(const int L, const int R);
inline int update(const int L);
inline void allocRecvArray(Arrays *send, Arrays *recv, 
			swInt diagNum, swInt offDiagNum);
inline void outputArray(Arrays *slave, Arrays *host, swInt startIdx);

void directSegmentIterator_e2v_slave(DS_edge2VertexPara *para)
{
// **********************************************************************
// Start Initiation
// TODO
// move to another function\
// **********************************************************************
	myId = athread_get_id(-1);

	swInt *segStarts,*edgeStarts,*owner,*neighbor;
	segStarts  = para->segStarts;
	edgeStarts = para->edgeStarts;
	owner      = para->startVertex;
	neighbor   = para->endVertex;

	swFloat *upper,*lower,*diag,*vertex;
	upper  = para->frontEdgeData;
	lower  = para->backEdgeData;
	diag   = para->selfConnData;
	vertex = para->vertexData;

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
	swInt colNum[BLOCKNUM64K+1],accuColNum[BLOCKNUM64K+2];
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
	swInt colIdx[edgeLen],edgeGlobalIdx[edgeLen];
//	swInt neighbor_slave[maxEdges],sNeighbor_slave[maxEdges];
	owner_slave     = (swInt*)ldm_malloc(maxEdges*sizeof(swInt));
	neighbor_slave  = (swInt*)ldm_malloc(maxEdges*sizeof(swInt));
	sNeighbor_slave = (swInt*)ldm_malloc(maxEdges*sizeof(swInt));
    DMA_Get(&sNeighbor_slave[0],&neighbor[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
    DMA_Get(&owner_slave[0],&owner[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));

    for(i=0;i<BLOCKNUM64K+1;i++) {colNum[i]=0;}
	for(i=0;i<edgeLen;i++)
	{
		if(sNeighbor_slave[i]<segStarts_slave[BLOCKNUM64K])
		{
			for(j=myId;j<BLOCKNUM64K;j++)
			{
				if(sNeighbor_slave[i]<segStarts_slave[j+1]
							&& sNeighbor_slave[i]>=segStarts_slave[j])
				{
					colIdx[i] = j;
					colNum[j]++;
				}
			}
		} else
		{
			colIdx[i]=BLOCKNUM64K;
			colNum[BLOCKNUM64K]++;
		}
	}

    accuColNum[0]=0;
	for(i=0;i<BLOCKNUM64K+1;i++)
	{
		accuColNum[i+1] += accuColNum[i]+colNum[i];
	}
#ifdef DEBUG
	assert(accuColNum[myId+1]==colNum[myId]);
	assert(accuColNum[BLOCKNUM64K+1]==edgeLen);
#endif

    // Get the global index of each edge for quicksort
	// priority: column block index > row index > col index
    for(i=0;i<edgeLen;i++)
	{
		edgeGlobalIdx[i]
			= colIdx[i]*cellLen*vertexNumber
			+ owner_slave[i]*vertexNumber
			+ sNeighbor_slave[i];
}
// ***********************************************************************
// END Initiation
// ***********************************************************************

	// data stored in slave core
//	sVertex_slave   = (Arrays*)ldm_malloc(sizeof(Arrays));
//	lowVertex_slave = (Arrays*)ldm_malloc(sizeof(Arrays));
//	upVertex_slave  = (Arrays*)ldm_malloc(sizeof(Arrays));
//	upper_slave     = (Arrays*)ldm_malloc(sizeof(Arrays));
//	sLower_slave    = (Arrays*)ldm_malloc(sizeof(Arrays));
//	rLower_slave    = (Arrays*)ldm_malloc(sizeof(Arrays));
//	diag_slave      = (Arrays*)ldm_malloc(sizeof(Arrays));
//
//	Arrays sVertex_slave,lowVertex_slave,upVertex_slave,upper_slave;
//	Arrays sLower_slave,rLower_slave,diag_slave;
//	if(myId==0) printf("vertex-fArrayNum: %d\n",vertex.fArrayNum);
    copyArray(&sVertex_slave, vertex, cellLen, segStarts_slave[myId]);
	copyArray(&upper_slave,   upper,  edgeLen, edgeStarts_slave[myId]);
	copyArray(&sLower_slave,  lower,  edgeLen, edgeStarts_slave[myId]);
	copyArray(&diag_slave,    diag,   cellLen, segStarts_slave[myId]);

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
    pivLower
		= (swFloat**)ldm_malloc(sLower_slave.fArrayNum*sizeof(swFloat*));
	pivUpper
		= (swFloat**)ldm_malloc(upper_slave.fArrayNum*sizeof(swFloat*));
	for(i=0;i<sLower_slave.fArrayNum;i++)
	{
		pivLower[i] = (swFloat*)ldm_malloc
			(sLower_slave.fArrayDims[i]*sizeof(swFloat));
	}
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		pivUpper[i] = (swFloat*)ldm_malloc
			(upper_slave.fArrayDims[i]*sizeof(swFloat));
	}

//if(myId==0) printArray("%d",sNeighbor_slave,edgeLen);
	quickSort(&edgeGlobalIdx[0],edgeLen);

	// Communicate neighbor
    swInt ownNeiSendIdx[BLOCKNUM64K],ownNeiRecvIdx[BLOCKNUM64K];
	swInt ownNeiRecvList[BLOCKNUM64K];
	int lastId,startIdx;
    for(i=0;i<colNum[myId];i++)
	{
		neighbor_slave[i] = sNeighbor_slave[i];
	}
	for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
	for(i=0;i<_total_send_pcg;i++)
	{
		if(myId>_sPacks[i].dst_id) continue;
		startIdx
			=accuColNum[_sPacks[i].dst_id]+ownNeiSendIdx[_sPacks[i].dst_id];
		for(j=0;j<_sPacks[i].cva;j++)
		{
			_sPacks[i].data[j] = sNeighbor_slave[startIdx+j];
		}
		ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].cva;
	}
	transform_data();
	recvNum=colNum[myId];
	lastId=0;
	int k=0;
	for(i=0;i<BLOCKNUM64K;i++) {ownNeiRecvList[i]=colNum[myId];}
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
	recvNum-=colNum[myId];
//if(myId==5) printArray("%d",neighbor_slave,colNum[myId]+recvNum);
//if(myId==0) printArray("%d",colIdx,edgeLen);
#ifdef DEBUG
	assert(maxEdges>recvNum+colNum[myId]);
#endif

	allocRecvArray(&rLower_slave, &sLower_slave,  colNum[myId], recvNum);
	allocRecvArray(&lowVertex_slave, &sVertex_slave, cellLen, recvNum);
	allocRecvArray(&upVertex_slave, &sVertex_slave,
				cellLen, edgeLen-colNum[myId]);

	// Communicate Lower
	int iArray,dims,iDim;
    for(iArray=0;iArray<rLower_slave.fArrayNum;iArray++)
	{
		dims = rLower_slave.fArrayDims[iArray];
		for(iDim=0;iDim<dims;iDim++)
		{
			for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
			for(i=0;i<_total_send_pcg;i++)
			{
				if(myId>_sPacks[i].dst_id) continue;
				startIdx = accuColNum[_sPacks[i].dst_id] * dims
					+ ownNeiSendIdx[_sPacks[i].dst_id] + iDim;
				for(j=0;j<_sPacks[i].cva;j++)
				{
					_sPacks[i].data[j]
						=sLower_slave.floatArrays[iArray][startIdx+j*dims];
				}
				ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].cva*dims;
			}
			transform_data();
			startIdx=colNum[myId]*dims+iDim;
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
						= accuColNum[_sPacks[i].dst_id]
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

	for(i=0;i<maxEdges;i++)
	{
		owner_slave[i]     -= segStarts_slave[myId];
		sNeighbor_slave[i] -= segStarts_slave[myId];
		neighbor_slave[i]  -= segStarts_slave[myId];
	}
	for(i=colNum[myId];i<maxEdges;i++)
	{
		sNeighbor_slave[i] = cellLen+i-colNum[myId];
	}
//if(myId==1) printArray("%d",sNeighbor_slave,edgeLen);
//if(myId==1) printf("cellLen: %d, colNum[myId]: %d\n",cellLen,colNum[myId]);

	// Compute the upper triangle
    rLower_slave.fArraySizes    = 0;
	diag_slave.fArraySizes  = 0;
//	para->fun_slave(&upper_slave,&rLower_slave,&diag_slave,&upVertex_slave,
//				owner_slave,sNeighbor_slave);
	//selfConn computation
	diag	= accessArray(&diag_slave, 0);
	swFloat* x		= accessArray(&upVertex_slave, 0);
	swFloat* b		= accessArray(&upVertex_slave, 1);
	//frontEdge computation
	upper	= accessArray(&upper_slave, 0);
	swInt edgeNumber = getArraySize(&upper_slave);
	swInt iedge;
	for( iedge = 0; iedge < edgeNumber; iedge++)
{
//if(owner_slave[iedge]+segStarts_slave[myId]==73) printf("upper:%d,%d,%d,%f,%f,%f\n",iedge,owner_slave[iedge],sNeighbor_slave[iedge],b[owner_slave[iedge]],upper[iedge],x[sNeighbor_slave[iedge]]);
			b[owner_slave[iedge]] += upper[iedge]*x[sNeighbor_slave[iedge]];
}

//if(myId==4) printArray("%f",upVertex_slave.floatArrays[1],cellLen);

	// Compute the lower triangle
	for(i=colNum[myId];i<maxEdges;i++)
	{
		owner_slave[i] = cellLen+i-colNum[myId];
	}
	x		= accessArray(&lowVertex_slave, 0);
	b		= accessArray(&lowVertex_slave, 1);
    upper_slave.fArraySizes  = 0;
	rLower_slave.fArraySizes = recvNum+colNum[myId];
//	para->fun_slave(&upper_slave,&rLower_slave,&diag_slave,
//				&lowVertex_slave,owner_slave,neighbor_slave);
	//backEdge computation
	lower	= accessArray(&rLower_slave, 0);
	edgeNumber = getArraySize(&rLower_slave);
	for( iedge = 0; iedge < edgeNumber; iedge++)
	{
//if(neighbor_slave[iedge]+segStarts_slave[myId]==73) printf("lower:%d,%d,%d,%f,%f,%f\n",iedge,owner_slave[iedge],neighbor_slave[iedge],b[neighbor_slave[iedge]],lower[iedge],x[owner_slave[iedge]]);
			b[neighbor_slave[iedge]] += lower[iedge]*x[owner_slave[iedge]];
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
//if(ivertex+segStarts_slave[myId]==73) printf("diag:%f,%f,%f\n",b[ivertex],diag[ivertex],x[ivertex]);
		b[ivertex] += diag[ivertex]*x[ivertex];
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
	outputArray(&diag_slave,      diag,   segStarts_slave[myId]);
	outputArray(&lowVertex_slave, vertex, segStarts_slave[myId]);
	outputArray(&upper_slave,     upper,  edgeStarts_slave[myId]);
}


inline void copyArray(Arrays *slave, Arrays *host,
			            swInt length, swInt startIdx)
{
	int i;
//	if(myId==0) printf("%d,%d\n",length,host->fArrayNum);
	swInt fArrayNum    = host->fArrayNum;
	slave->fArrayNum   = fArrayNum;
	slave->fArraySizes = length;
	slave->fArrayDims  = (swInt)ldm_malloc(fArrayNum*sizeof(swInt));
	slave->fArrayInOut = (swInt)ldm_malloc(fArrayNum*sizeof(swInt));
	DMA_Get(slave->fArrayDims,host->fArrayDims,fArrayNum*sizeof(swInt));
	DMA_Get(slave->fArrayInOut,host->fArrayInOut,fArrayNum*sizeof(swInt));
	slave->floatArrays = (swFloat**)ldm_malloc(fArrayNum*sizeof(swFloat*));
	for(i=0;i<fArrayNum;i++)
	{
		slave->floatArrays[i] = (swFloat*)ldm_malloc
			(length*slave->fArrayDims[i]*sizeof(swFloat));
		DMA_Get(slave->floatArrays[i],
			&host->floatArrays[i][slave->fArrayDims[i]*startIdx],
			length*slave->fArrayDims[i]*sizeof(swFloat));
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
	pivOwner    = owner_slave[L];
	pivNeighbor = sNeighbor_slave[L];
	for(i=0;i<sLower_slave.fArrayNum;i++)
	{
		for(j=0;j<sLower_slave.fArrayDims[i];j++)
		  pivLower[i][j] = sLower_slave.floatArrays[i][L+j];
	}
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		for(j=0;j<upper_slave.fArrayDims[i];j++)
		  pivUpper[i][j] = upper_slave.floatArrays[i][L+j];
	}
}

inline int exchange(const int L, const int R)
{
	int i,j;
	owner_slave[L]     = owner_slave[R];
	sNeighbor_slave[L] = sNeighbor_slave[R];
	for(i=0;i<sLower_slave.fArrayNum;i++)
	{
		for(j=0;j<sLower_slave.fArrayDims[i];j++)
		  sLower_slave.floatArrays[i][L+j]
			  = sLower_slave.floatArrays[i][R+j];
	}
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		for(j=0;j<upper_slave.fArrayDims[i];j++)
		  upper_slave.floatArrays[i][L+j]
			  = upper_slave.floatArrays[i][R+j];
	}
}

inline int update(const int L)
{
	int i,j;
	owner_slave[L]     = pivOwner;
	sNeighbor_slave[L] = pivNeighbor;
	for(i=0;i<sLower_slave.fArrayNum;i++)
	{
		for(j=0;j<sLower_slave.fArrayDims[i];j++)
		  sLower_slave.floatArrays[i][L+j] = pivLower[i][j];
	}
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		for(j=0;j<upper_slave.fArrayDims[i];j++)
		  upper_slave.floatArrays[i][L+j] = pivUpper[i][j];
	}
}

inline void allocRecvArray(Arrays *recv, Arrays *send, 
			swInt diagNum, swInt offDiagNum)
{
	int i,j;
	recv->fArrayNum   = send->fArrayNum;
	recv->fArraySizes = diagNum+offDiagNum;
	recv->fArrayDims  = (swInt)ldm_malloc(recv->fArrayNum*sizeof(swInt));
	recv->fArrayInOut = (swInt)ldm_malloc(recv->fArrayNum*sizeof(swInt));
	recv->floatArrays
		= (swFloat**)ldm_malloc(recv->fArrayNum*sizeof(swFloat*));
	for(i=0;i<recv->fArrayNum;i++)
	{
		recv->fArrayDims[i] = send->fArrayDims[i];
		recv->fArrayInOut[i] = send->fArrayInOut[i];
		recv->floatArrays[i]
			=(swFloat*)ldm_malloc
			((diagNum+offDiagNum)*recv->fArrayDims[i]*sizeof(swFloat));
		for(j=0;j<diagNum*recv->fArrayDims[i];j++)
		{
			recv->floatArrays[i][j]
				= send->floatArrays[i][j];
		}
	}
}

inline void outputArray(Arrays *slave, Arrays *host, swInt startIdx)
{
	int i;
	swInt length = slave->fArraySizes;
	for(i=0;i<slave->fArrayNum;i++)
	{
		if(slave->fArrayInOut[i]==COPYIN) continue;
		DMA_Put(&(host->floatArrays[i][startIdx*slave->fArrayDims[i]]),
					slave->floatArrays[i],
					length*slave->fArrayDims[i]*sizeof(swFloat));
	}
}


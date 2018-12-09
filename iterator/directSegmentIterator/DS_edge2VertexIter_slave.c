#include <stdlib.h>
#include <stdio.h>
#include "slave.h"
#include "simd.h"
#include "dma.h"
#include "rlmpi.h"
#include "swMacro.h"
#include "iterator_struct.h"
#include "userFunc_slave.h"
#include "DS_edge2VertexIter_host.h"

__thread_local volatile int myId;

//extern void (*operatorFunPointer_s)(MLBFunParameters *MLBFunParas);
//extern swInt spIndex, maxXNum, maxCell, maxEdge, cpeBlockNum;
//extern swInt *blockStarts, *ownNeiSendList, *owner, *neighbor;
//extern swFloat *upper, *lower, *diag, *x, *b;
//extern int isXExist;

void DS_edge2VertexIter_slave(DS_e2VParameters *e2VParas)
{
	segStarts   = e2VParas->DSParas->segStarts;
	edgeStarts  = e2VParas->DSParas->edgeStarts;
	owner       = e2VParas->DSParas->owner;
	neighbor    = e2VParas->DSParas->neighbor;

	upper  = e2VParas->upper;
	lower  = e2VParas->lower;
	diag   = e2VParas->diag;
	vertex = e2VParas->vertex;

	spIndex      = e2VParas->spIndex;
	isXExist     = e2VParas->isXExist;
	maxEdges     = e2VParas->DSParas->maxEdges;
	vertexNumber = e2VParas->DSParas->vertexNumber;
	edgeNumber   = e2VParas->DSParas->edgeNumber;
	
	swInt i,j,k,row,col,maxEdges,startIdx,lastId,edgeLen,cellLen,recvNum;
	myId = athread_get_id(-1);
	row = myId + spIndex*BLOCKNUM64K;
	swInt ownNeiSendIdx[BLOCKNUM64K],ownNeiRecvIdx[BLOCKNUM64K];
	swInt ownNeiRecvList[BLOCKNUM64K];

//	volatile swInt total_send_pcg = _table_ldm.nPUTR + _table_ldm.nPutrSameRow + _table_ldm.nPutcSameCol;
//	volatile swInt totalLength
//		= maxXNum*3*sizeof(swFloat)
//		+ maxEdge*2*(sizeof(swFloat)+sizeof(swInt))
//		+ maxCell*2*sizeof(swInt)
//		+ maxCell*sizeof(swFloat)
//		+ BLOCKNUM64K*8*sizeof(swInt)
//		+ total_send_pcg*2*32;
//	if(totalLength > 64*1024*0.95 && myId==0)
//	{
////	if(myId == 0){
//		printf("The allocated LDM exceeds 64KB, the memory size is %d\n",
//					totalLength);
//	}

	segStarts_slave  = (swInt*)ldm_malloc((BLOCKNUM64K+1)*sizeof(swInt));
	DMA_Get(&segStarts_slave[0],&segStarts[spIndex*BLOCKNUM64K],
				(BLOCKNUM64K+1)*sizeof(swInt));
	edgeStarts_slave = (swInt*)ldm_malloc((BLOCKNUM64K+1)*sizeof(swInt));
	DMA_Get(&edgeStarts_slave[0],&edgeStarts[spIndex*BLOCKNUM64K],
				(BLOCKNUM64K+1)*sizeof(swInt));
	colNum     = (swInt*)ldm_malloc((BLOCKNUM64K+1)*sizeof(swInt));
	accuColNum = (swInt*)ldm_malloc((BLOCKNUM64K+2)*sizeof(swInt));
	cellLen    = segStarts_slave[myId+1]-segStarts_slave[myId];

	// Get the maximum edge count of all segments.
	// assert(maxEdges>recvNum+sendNum);
	maxEdges = 0;
	for(i=0;i<BLOCKNUM64K;i++)
	{
		edgeLen  = edgeStarts_slave[i+1]-edgeStarts_slave[i];
		maxEdges = maxEdges > edgeLen ? maxEdges : edgeLen;
	}
	maxEdges = maxEdges*1.02;

//	x_slave = (swloat*)ldm_malloc(cellLen*sizeof(swFloat));
//	b_slave = (swloat*)ldm_malloc(cellLen*sizeof(swFloat));

	edgeLen        = edgeStarts_slave[myId+1]-edgeStarts_slave[myId];
	colIdx         = (swInt*)ldm_malloc(edgeLen*sizeof(swInt));
	colGlobalIdx   = (swInt*)ldm_malloc(edgeLen*sizeof(swInt));
	owner_slave    = (swInt*)ldm_malloc(maxEdges*sizeof(swInt));
	neighbor_slave = (swInt*)ldm_malloc(maxEdges*sizeof(swInt));

	sVertex_slave = (Arrays*)ldm_malloc(sizeof(Arrays));
	rVertex_slave = (Arrays*)ldm_malloc(sizeof(Arrays));
	upper_slave   = (Arrays*)ldm_malloc(sizeof(Arrays));
	sLower_slave  = (Arrays*)ldm_malloc(sizeof(Arrays));
	rLower_slave  = (Arrays*)ldm_malloc(sizeof(Arrays));
	diag_slave    = (Arrays*)ldm_malloc(sizeof(Arrays));
	copyArrays(sVertex_slave, vertex, cellLen, segStarts_slave[myId]);
	copyArrays(upper_slave,   upper,  edgeLen, edgeStarts_slave[myId]);
	copyArrays(sLower_slave,  lower,  edgeLen, edgeStarts_slave[myId]);
	copyArrays(diag_slave,    diag,   cellLen, segStarts_slave[myId]);


//	DMA_Get(&b_slave[0],&b[segStarts_slave[myId]],cellLen*sizeof(swFloat));
//	DMA_Get(&x_slave[0],&x[segStarts_slave[myId]],cellLen*sizeof(swFloat));
	swInt *sNeighbor_slave,*rNeighbor_slave,*sOwner_slave;
	sOwner_slave    = (swInt*)ldm_malloc(edgeLen*sizeof(swInt));
	sNeighbor_slave = (swInt*)ldm_malloc(edgeLen*sizeof(swInt));
	rNeighbor_slave = (swInt*)ldm_malloc((maxEdges-edgeLen)*sizeof(swInt));
	DMA_Get(&sOwner_slave[0],*owner[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
	DMA_Get(&sNeighbor_slave[0],*neighbor[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
//	DMA_Get(&sLower_slave[0],*lower[edgeStarts_slave[myId]],
//				edgeLen*sizeof(swFloat));
//	DMA_Get(&sUpper_slave[0],*upper[edgeStarts_slave[myId]],
//				edgeLen*sizeof(swFloat));

	FunParameters FunParas = {sLower_slave,sUpper_slave,diag_slave,
		vertex_slave,sOwner_slave,sNeighbor_slave,0,0,0,0};

	// Compute the diagonal data
	// Get the column block index of each edge: colIdx
	// Get the count of the edges in each column block: colNum
	for(i=0;i<BLOCKNUM64K+1;i++) {colNum[i]=0;}
	for(i=0;i<edgeLen;i++)
	{
		if(sNeighbor_slave[i]<segStarts_slave[myId+1]
					&& sNeighbor_slave[i]>=segStarts_slave[myId])
		{
			// 计算对角块非对角元素
			// 调用函数指针：flag=0
//			b_slave[sOwner_slave[i]]
//				+= sUpper_slave[i]*x_slave[sNeighbor_slave[i]];
//			b_slave[sNeighbor_slave[i]]
//				+= sLower_slave[i]*x_slave[sOwn_slave[i]];
			colIdx[i] = BLOCKNUM64K-myId;
			cellNum[myId]++;
		} else if(sNeighbor_slave[i]<segStarts_slave[BLOCKNUM64K])
		{
			for(j=myId+1;j<BLOCKNUM64K;j++)
			{
				if(sNeighbor_slave[i]<segStarts_slave[j+1]
						&& sNeighbor_slave[i]>=segStarts_slave[j])
				{
					colIdx[i] = BLOCKNUM64K-j;
					colNum[j]++;
				}
			}
		} else
		{
			// TODO
			// The start/end neighbor index need to be recorded 
			// if the off-diagonal block is computed in slave core;
			colIdx[i]=0;
			colNum[BLOCKNUM64K]++;
		}
	}
	accuColNum[0]=0;
	for(i=0;i<BLOCKNUM64K+1;i++)
	{
		accuColNum[i+1] += accuColNum[i]+colNum[i];
	}
	// Get the global index of each edge for quicksort
	// priority: column block index > row index > col index
	for(i=0;i<edgeLen;i++)
	{
		edgeGlobalIdx[i]
			= colIdx[i]*cellLen*vertexNumber
			+ sOwner_slave[i]*vertexNumber
			+ sNeighbor_slave[i];
	}

	// Record the compared base in quicksort.
	pivLower
		= (swFloat**)ldm_malloc(sLower_slave->fArrayNum*sizeof(swFloat*));
	pivUpper
		= (swFloat**)ldm_malloc(sUpper_slave->fArrayNum*sizeof(swFloat*));
	for(i=0;i<sLower_slave->fArrayNum;i++)
	{
		pivLower[i]
			= (swFloat*)ldm_malloc(sLower->fArrayDims[i]*sizeof(swFloat));
	}
	for(i=0;i<sUpper_slave->fArrayNum;i++)
	{
		pivUpper[i]
			= (swFloat*)ldm_malloc(sUpper->fArrayDims[i]*sizeof(swFloat));
	}

	quickSort(&edgeGlobalIdx[0],edgeLen);

	// Communicate neighbor
	for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
	for(i=0;i<_total_send_pcg;i++)
	{
		if(myId>_sPacks[i].dst_id) continue;
		startIdx
		   = accuColNum[_sPacks[i].dst_id]+ownNeiSendIdx[_sPacks[i].dst_id];
		for(j=0;j<_sPacks[i].indM;j++)
		{
			_sPacks[i].data[j] = sNeighbor_slave[startIdx+j];
		}
		ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].indM;
	}
	transform_data();
	recvNum=0;
	lastId=0;
	for(i=0;i<BLOCKNUM64K;i++) {ownNeiRecvIdx[i]=0;}
	for(i=0;i<_total_send_pcg;i++)
	{
		if(_rPacks[i].src_id>myId) continue;
		for(j=0;j<_rPacks[i].indM;j++)
		{
			rNeighbor_slave[recvNum] = _rPacks[i].data[j];
			recvNum++;
		}
		if(_rPacks[i].src_id!=lastId) {k++;lastId=_rPacks[i].src_id;}
		ownNeiRecvList[k+1]=lowerNum;
		ownNeiRecvIdx[_rPacks[i].src_id]=k;
	}

	// reallocate the lower data and vertex data for register communication
	allocRecvArrays(edgeLen-colNum[myId], recvNum);

	// Communicate Lower
	for(iArray=0;iArray<rLower_slave->fArrayNum;iArray++)
	{
		for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
		for(i=0;i<_total_send_pcg;i++)
		{
			if(myId>_sPacks[i].dst_id) continue;
			startIdx
				= accuColNum[_sPacks[i].dst_id]
				+ ownNeiSendIdx[_sPacks[i].dst_id];
			for(j=0;j<_sPacks[i].indM;j++)
			{
				_sPacks[i].data[j]
					= sLower_slave[iArray][startIdx+j];
//				* x_slave[sOwner_slave[startIdx+j]-segStarts_slave[myId]];
			}
			ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].indM;
		}
		transform_data();
		startIdx=0;
		for(i=0;i<_total_send_pcg;i++)
		{
			if(_rPacks[i].src_id>myId) continue;
			for(j=0;j<_rPacks[i].indM;j++)
			{
				rLower_slave[iArrays][startIdx] = _rPacks[i].data[j];
				startIdx++;
			}
		}
	}

//	for(i=0;i<lower_slave->fArrayNum;i++)
//	{
//		for(j=0;j<sendNum;j++)
//		{
//			lower_slave->floatArrays[i][recvNum+j]
//				= upper_slave->floatArrays[i][j];
//		}
//	}
	// Compute the lower triangle
//	for(i=0;i<startIdx;i++)
//	{
//		b_slave[rNeighbor_slave[i]-segStarts_slave[myId]]+=rLower_slave[i];
//	}

	// Communicate X
	for(iArray=0;iArray<sVertex_slave->fArrayNum;iArray++)
	{
		if(sVertex_slave->fArrayInOut[iArray]==COPYOUT) continue;
		for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
		for(i=0;i<_total_send_pcg;i++)
		{
			if(myId<_sPacks[i].dst_id)
			{
				startIdx
					= accuColNum[_sPacks[i].dst_id]
					+ ownNeiSendIdx[_sPacks[i].dst_id];
				for(j=0;j<_sPacks[i].indM;j++)
				{
					_sPacks[i].data[j]
						= sVertex_slave[iArray]
						[sOwner_slave[startIdx+j]-segStarts_slave[myId]];
				}
			} else
			{
				startIdx
				   = ownNeiRecvList[ownNeiRecvIdx[_sPacks[i].dst_id]]
				   + ownNeiSendIdx[_sPacks[i].dst_id];
				for(j=0;j<_sPacks[i].indM;j++)
				{
					_sPacks[i].data[j]
						= sVertex_slave[iArray]
						[rNeighbor_slave[startIdx+j]-segStarts_slave[myId]];
				}
			}
			ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].indM;
		}
		transform_data();
		startIdx=0;
		for(i=0;i<_total_send_pcg;i++)
		{
	//		if(_rPacks[i].src_id>myId) continue;
			for(j=0;j<_rPacks[i].indM;j++)
			{
				rVertex_slave[iArray][startIdx] = _rPacks[i].data[j];
				startIdx++;
			}
		}
	}

	swInt nonZeroNum = recvNum + edgeLen - colNum[myId];
	for(i=0;i<recvNum;i++)
	{
		neighbor_slave[i] = i;
		owner_slave[i]    = rNeighbor_slave[i];
	}
	for(i=0;i<edgeLen-colNum[myId];i++)
	{
		neighbor_slave[i] = i;
		owner_slave[i]    = sOwner_slave[i+colNum[myId]];
	}

	FunParas->owner    = &owner_slave[0];
	FunParas->neighbor = &neighbor_slave[0];
	FunParas->lower    = rLower_slave;
	FunParas->vertex   = rVertex_slave;
	// Compute the off-diagonal block
//	for(i=0;i<startIdx;i++)
//	{
//		b_slave[sNeighbor_slave[i]-segStarts_slave[myId]]
//			= sLower_slave[i]*recvX_slave[i];
//	}

	// TODO
	// Assume that lower data won't be output
//	outputArray(rLower_slave);
	outputArray(diag_slave);
	outputArray(rVertex_slave);
	outputArray(upper_slave);
//	DMA_Put(&b[segStarts_slave[myId]],&b_slave[0],cellLen*sizeof(swFloat));
}

inline void copyArray(Arrays *slave, Arrays *host,
			swInt length, swInt startIdx)
{
	swInt fArrayNum = host->fArrayNum;
	slave->fArrayNum = fArrayNum;
	slave->fArraySize = length;
	DMA_Get(slave->fArrayDims,host->fArrayDims,fArrayNum*sizeof(swInt));
//	DMA_Get(slave->fArraySizes,host->fArraySizes,fArrayNum*sizeof(swInt));
	DMA_Get(slave->fArrayInOut,host->fArrayInOut,fArrayNum*sizeof(swInt));
	slave->floatArrays = (swFloat**)ldm_malloc(fArrayNum*sizeof(swFloat*));
	for(i=0;i<fArrayNum;i++)
	{
		slave->floatArrays[i]
			= (swFloat*)ldm_malloc
			(length*slave->fArrayDims[i]*sizeof(swFloat));
		DMA_Get(slave->floatArrays[i],
			&host->floatArrays[i][startIdx],
			length*slave->fArrayDims[i]*sizeof(swFloat));
	}
}

inline void allocRecvArray(swInt sendNum, swInt recvNum)
{
	int i,j;
	rLower_slave->fArrayNum  = sLower_slave->fArrayNum;
	rLower_slave->fArraySize = recvNum;
	for(i=0;i<rLower_slave->fArrayNum;i++)
	{
		rLower_slave->fArrayDims[i] = sLower_slave->fArrayDims[i];
		rLower_slave->fArrayInOut[i] = sLower_slave->fArrayInOut[i];
		rLower_slave->floatArrays[i]
			=(swFloat*)ldm_malloc(recvNum*sizeof(swFloat));
	}

	swInt cellLen = segStarts_slave[myId+1]-segStarts_slave[myId];
	rVertex_slave->fArrayNum = sVertex_slave->fArrayNum;
	rVertex_slave->fArraySize = sendNum+recvNum;
	for(i=0;i<rVertex_slave->fArrayNum;i++)
	{
		rVertex_slave->fArrayDims[i] = sVertex_slave->fArrayDims[i];
		rVertex_slave->fArrayInOut[i] = sVertex_slave->fArrayInOut[i];
		rVertex_slave->floatArrays[i]
			=(swFloat*)ldm_malloc((sendNum+recvNum)*sizeof(swFloat));
		if(rVertex_slave->fArrayInOut[i]==COPYOUT)
		{
			for(j=0;j<cellLen;j++)
			{
				rVertex_slave->floatArrays[i][j]
					= sVertex_slave->floatArrays[i][j];
			}
		}
	}
}

inline void outputArray(Arrays *slave, Arrays *host)
{
	int i;
	swInt length = slave->fArraySize;
	for(i=0;i<slave->fArrayNum;i++)
	{
		if(slave->fArrayInOut[i]=COPYIN) continue;
		DMA_Put(&host->floatArray[i][segStarts_slave[myId]],
				slave->floatArray[i],
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
				if (L<R) arr[L]=arr[R]; exange(L,R); L++;
				while (arr[L]<=piv && L<R) L++; 
				if (L<R) arr[R]=arr[L]; exange(R,L); R--;
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
	pivOwner    = sOwner_slave[L];
	pivNeighbor = sNeighbor_slave[L];
	for(i=0;i<sLower_slave->fArrayNum;i++)
	{
		for(j=0;j<sLower_slave->fArrayDims[i];j++)
			pivLower[i][j] = sLower_slave->floatArrays[i][L+j];
	}
	for(i=0;i<upper_slave->fArrayNum;i++)
	{
		for(j=0;j<upper_slave->fArrayDims[i];j++)
			pivUpper[i][j] = upper_slave->floatArrays[i][L+j];
	}
}

inline int exange(const int L, const int R)
{
	int i,j;
	sOwner_slave[L]    = sOwner_slave[R];
	sNeighbor_slave[L] = sNeighbor_slave[R];
	for(i=0;i<sLower_slave->fArrayNum;i++)
	{
		for(j=0;j<sLower_slave->fArrayDims[i];j++)
			sLower_slave->floatArrays[i][L+j]
				= sLower_slave->floatArrays[i][R+j];
	}
	for(i=0;i<upper_slave->fArrayNum;i++)
	{
		for(j=0;j<upper_slave->fArrayDims[i];j++)
			upper_slave->floatArrays[i][L+j]
				= upper_slave->floatArrays[i][R+j];
	}
}

inline int update(const int L)
{
	int i,j
	sOwner_slave[L]    = pivOwner;
	sNeighbor_slave[L] = pivNeighbor;
	for(i=0;i<sLower_slave->fArrayNum;i++)
	{
		for(j=0;j<sLower_slave->fArrayDims[i];j++)
			sLower_slave->floatArrays[i][L+j] = pivLower[i][j];
	}
	for(i=0;i<upper_slave->fArrayNum;i++)
	{
		for(j=0;j<upper_slave->fArrayDims[i];j++)
			upper_slave->floatArrays[i][L+j] = pivUpper[i][j];
	}
}

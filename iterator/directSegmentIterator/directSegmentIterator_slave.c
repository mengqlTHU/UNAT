#include "swMacro.h"
#include "slave.h"
#include <assert.h>
#include "rlmpi.h"
#include "slaveUtils.h"
#include "directSegmentIterator.h"
#include "iterator.h"
#include "onChipTrans.h"

#define DEBUGID 0
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

#define copyArray(slave, host, edgeLen, startIdx, recvEdges, DMA) \
{ \
	int i; \
	swInt fArrayNum    = host->fArrayNum; \
	slave.fArrayNum   = fArrayNum; \
	slave.fArraySizes = edgeLen+recvEdges; \
	if(fArrayNum>0) \
	{ \
	LDM_NEW(slave.fArrayDims, swInt, fArrayNum); \
	LDM_NEW(slave.fArrayInOut, swInt, fArrayNum); \
	DMA_Get(slave.fArrayDims,host->fArrayDims,fArrayNum*sizeof(swInt)); \
	DMA_Get(slave.fArrayInOut,host->fArrayInOut,fArrayNum*sizeof(swInt));\
	LDM_NEW(slave.floatArrays, swFloat*, fArrayNum); \
	for(i=0;i<fArrayNum;i++) \
	{ \
		LDM_NEW(slave.floatArrays[i], swFloat, \
					(edgeLen+recvEdges)*slave.fArrayDims[i]); \
		if(!DMA) continue; \
		if(edgeLen<=0) continue; \
		DMA_Get(&slave.floatArrays[i][recvEdges*slave.fArrayDims[i]], \
			&host->floatArrays[i][slave.fArrayDims[i]*startIdx], \
			edgeLen*slave.fArrayDims[i]*sizeof(swFloat)); \
	} \
	} \
}

#define initRegisterPacks(rlmpi_info) \
{  \
    _table_ldm.nGetcSkew = rlmpi_info->table[_MYID].nGetcSkew;\
    _table_ldm.nPutrSkew = rlmpi_info->table[_MYID].nPutrSkew;\
    _table_ldm.nGetrPutcSkew = rlmpi_info->table[_MYID].nGetrPutcSkew;\
    _nCycleSkew = rlmpi_info->nCycle;\
    int length = 32 * (sizeof (int8LDM) * _nCycleSkew / 32 + 1);\
    if (_nCycleSkew > 0) { \
		_putr_schedules_skew = (int8LDM*) ldm_malloc(length);\
		_getrputc_schedules_skew = (int8LDM*) ldm_malloc(length);\
		_getc_schedules_skew = (int8LDM*) ldm_malloc(length);\
		_get_reply = 0;\
		athread_get(PE_MODE, rlmpi_info->putr_schedules[_MYID], _putr_schedules_skew, length, &_get_reply, 0, 0, 0);\
		athread_get(PE_MODE, rlmpi_info->getrputc_schedules[_MYID], _getrputc_schedules_skew, length, &_get_reply, 0, 0, 0);\
		athread_get(PE_MODE, rlmpi_info->getc_schedules[_MYID], _getc_schedules_skew, length, &_get_reply, 0, 0, 0);\
		dma_wait(&_get_reply, 3);\
	} \
    _nCycleSameRow = rlmpi_info->nCycleSameRow;\
	_table_ldm.nGetrSameRow = rlmpi_info->table[_MYID].nGetrSameRow;\
	_table_ldm.nPutrSameRow = rlmpi_info->table[_MYID].nPutrSameRow;\
	length = 32 * (sizeof (int8LDM) * _nCycleSameRow / 32 + 1);\
	if (_nCycleSameRow > 0) { \
		_putr_schedules_same_row = (int8LDM*) ldm_malloc(length);\
		_getr_schedules_same_row = (int8LDM*) ldm_malloc(length);\
		_get_reply = 0;\
		athread_get(PE_MODE, rlmpi_info->putr_schedules_same_row [_MYID], _putr_schedules_same_row, length, &_get_reply, 0, 0, 0); \
		athread_get(PE_MODE, rlmpi_info->getr_schedules_same_row[_MYID], _getr_schedules_same_row, length, &_get_reply, 0, 0, 0); \
		dma_wait(&_get_reply, 2);\
	} \
    _nCycleSameCol = rlmpi_info->nCycleSameCol;\
	_table_ldm.nGetcSameCol = rlmpi_info->table[_MYID].nGetcSameCol;\
	_table_ldm.nPutcSameCol = rlmpi_info->table[_MYID].nPutcSameCol;\
	length = 32 * (sizeof (int8LDM) * _nCycleSameCol / 32 + 1);\
	if (_nCycleSameCol > 0) { \
		_putc_schedules_same_col = (int8LDM*) ldm_malloc(length);\
		_getc_schedules_same_col = (int8LDM*) ldm_malloc(length);\
		_get_reply = 0;\
		athread_get(PE_MODE, rlmpi_info->putc_schedules_same_col[_MYID], _putc_schedules_same_col, length, &_get_reply, 0, 0, 0);\
		athread_get(PE_MODE, rlmpi_info->getc_schedules_same_col[_MYID], _getc_schedules_same_col, length, &_get_reply, 0, 0, 0);\
		dma_wait(&_get_reply, 2);\
	} \
	load_rlmpi_data2(rlmpi_info);\
	ALLSYN;\
}

#define destroyRegisterPacks(rlmpi_info) \
{ \
	int length = 32 * (sizeof (int8LDM) * _nCycleSkew / 32 + 1);\
	if (_nCycleSkew > 0) {\
		ldm_free(_putr_schedules_skew, length);\
		ldm_free(_getrputc_schedules_skew, length);\
		ldm_free(_getc_schedules_skew, length);\
	}\
	_putr_schedules_skew = NULL;\
	_getrputc_schedules_skew = NULL;\
	_getc_schedules_skew = NULL;\
	\
    length = 32 * (sizeof (int8LDM) * _nCycleSameRow / 32 + 1);\
	\
	if (length > 0) {\
		ldm_free(_putr_schedules_same_row, length);\
		ldm_free(_getr_schedules_same_row, length);\
	}\
	length = 32 * (sizeof (int8LDM) * _nCycleSameCol / 32 + 1);\
	if (length > 0) {\
		ldm_free(_putc_schedules_same_col, length);\
		ldm_free(_getc_schedules_same_col, length);\
	}\
	if (_total_send_pcg > 0) {\
		ldm_free(_sPacks, _total_send_pcg * sizeof (Pack));\
	} \
	if (_total_recv_pcg > 0) {\
		ldm_free(_rPacks, _total_recv_pcg * sizeof (Pack));\
	}\
	ALLSYN;\
}

void directSegmentIterator_e2v_slave(DS_edge2VertexPara *para)
{
// **********************************************************************
// Start Initiation
// TODO
// move to another function\
// **********************************************************************

	DS_edge2VertexPara para_s;
	DMA_Get(&para_s,para,sizeof(DS_edge2VertexPara));
	swInt spIndex,vertexNumber;
	spIndex = para_s.spIndex;
	vertexNumber = para_s.vertexNumber;

	RlmpiInfo *rlmpi_info = &para_s.schedule_data[spIndex];
	initRegisterPacks(rlmpi_info);

	myId = athread_get_id(-1);

	int allocatable_ldm_size = 62000-_total_send_pcg*2*32;
//	if(myId==DEBUGID) printf("Allocatable LDM size: %d\n",allocatable_ldm_size);
	INIT_LDM_SPACE(44000);
	swInt *segStarts,*edgeStarts,*owner,*neighbor,*edgeNeiSeg,*accuColNum;
	swInt *rOwner,*rNeighbor,*colSegNum;
	segStarts  = para_s.segStarts;
	edgeStarts = para_s.edgeStarts;
	edgeNeiSeg = para_s.edgeNeiSeg;
	accuColNum = para_s.accuColNum;
	owner      = para_s.startVertex;
	neighbor   = para_s.endVertex;
	rOwner     = para_s.rOwner;
	rNeighbor  = para_s.rNeighbor;
	colSegNum  = para_s.colSegNum;

	Arrays *upper,*lower,*diag,*vertex,*rLower,*rUpper;
	upper  = para_s.frontEdgeData;
	lower  = para_s.backEdgeData;
	diag   = para_s.selfConnData;
	vertex = para_s.vertexData;
	rLower = para_s.rBackEdgeData;
	rUpper = para_s.rFrontEdgeData;


	Arrays  sVertex_slave,diag_slave;
//	Arrays  lowVertex_slave,upVertex_slave;
    swInt   *owner_slave,*neighbor_slave;
    swInt   *segStarts_slave;
    swInt   *accuColNum_slave,*edgeStarts_slave,*colSegNum_slave;
    swInt   *edgeNeiSeg_slave;
    Arrays  sLower_slave,upper_slave;
	swInt lowerFArrayNum = lower->fArrayNum;
	swInt diagFArrayNum  = diag->fArrayNum;
//if(myId==0) printf("%d,%d",lowerFArrayNum,diagFArrayNum);
	LDM_NEW(sLower_slave.floatArrays, swFloat*, lowerFArrayNum);
	LDM_NEW(sLower_slave.fArrayDims,  swFloat,  lowerFArrayNum);
	LDM_NEW(diag_slave.floatArrays,   swFloat*, diagFArrayNum);
	LDM_NEW(diag_slave.fArrayDims,    swFloat,  diagFArrayNum);

	volatile swInt cellLen,edgeLen;
	int i,j,row;
	row = spIndex*BLOCKNUM64K+myId;
	LDM_NEW(segStarts_slave,swInt,BLOCKNUM64K+1);
	LDM_NEW(edgeStarts_slave,swInt,BLOCKNUM64K+1);
	LDM_NEW(colSegNum_slave,swInt,BLOCKNUM64K);
	DMA_Get(&segStarts_slave[0],&segStarts[spIndex*BLOCKNUM64K],
				(BLOCKNUM64K+1)*sizeof(swInt));
	DMA_Get(&edgeStarts_slave[0],&edgeStarts[spIndex*BLOCKNUM64K],
				(BLOCKNUM64K+1)*sizeof(swInt));
	DMA_Get(colSegNum_slave,&colSegNum[spIndex*BLOCKNUM64K*BLOCKNUM64K+myId*BLOCKNUM64K],BLOCKNUM64K*sizeof(swInt));
//	LDM_NEW(colNum,swInt,BLOCKNUM64K+1);
	LDM_NEW(accuColNum_slave,swInt,BLOCKNUM64K+2);
	cellLen    = segStarts_slave[myId+1]-segStarts_slave[myId];
	DMA_Get(&accuColNum_slave[1],&accuColNum[spIndex*BLOCKNUM64K*(BLOCKNUM64K+1)+myId*(BLOCKNUM64K+1)],(BLOCKNUM64K+1)*sizeof(swInt));
    edgeLen        = edgeStarts_slave[myId+1]-edgeStarts_slave[myId];

    // Get the maximum edge count of all segments.
	// assert(maxEdges>recvNum+sendNum);
    swInt maxEdges,recvEdges;
//	for(i=0;i<BLOCKNUM64K;i++)
//	{
//		edgeLen  = edgeStarts_slave[i+1]-edgeStarts_slave[i];
//		maxEdges = maxEdges > edgeLen ? maxEdges : edgeLen;
//	}
//	maxEdges = maxEdges*1.5;
//	recvEdges = accuColNum_slave[myId+2]
//		+ para->recvStarts[myId+spIndex*BLOCKNUM64K];
	recvEdges = para_s.recvStarts[myId+spIndex*BLOCKNUM64K];

	int isLowerExisted = 1;
	if(lowerFArrayNum==0) isLowerExisted = 0;
//	maxEdges = MAX(recvEdges,edgeLen);

	LDM_NEW(owner_slave,swInt,edgeLen+recvEdges);
	LDM_NEW(neighbor_slave,swInt,edgeLen+recvEdges);
//	LDM_NEW(sNeighbor_slave,swInt,edgeLen);
	LDM_NEW(edgeNeiSeg_slave,swInt,edgeLen);
//    DMA_Get(&neighbor_slave[recvEdges],&neighbor[edgeStarts_slave[myId]],
//				edgeLen*sizeof(swInt));
////if(myId==DEBUGID) printArray("%d",neighbor_slave,edgeLen+recvEdges);
//    DMA_Get(&owner_slave[recvEdges],&owner[edgeStarts_slave[myId]],
//				edgeLen*sizeof(swInt));

	if(edgeLen>0) {
    DMA_Get(edgeNeiSeg_slave,&edgeNeiSeg[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
	}
//    DMA_Get(neighbor_slave,&neighbor[edgeStarts_slave[myId]],
//				edgeLen*sizeof(swInt));
//if(myId==DEBUGID) printArray("%d",neighbor_slave,edgeLen+recvEdges);

//    for(i=0;i<BLOCKNUM64K+1;i++) {colNum[i]=0;}
//if(myId==DEBUGID) printArray("%d",edgeNeiSeg_slave,edgeLen);
//	for(i=0;i<edgeLen;i++)
//	{
//		if(neighbor_slave[i+recvEdges]<segStarts_slave[BLOCKNUM64K])
//		{
//			if()
//			edgeNeiSeg_slave[i] -= spIndex*BLOCKNUM64K;
////			for(j=myId;j<BLOCKNUM64K;j++)
////			{
////				if(sNeighbor_slave[i]<segStarts_slave[j+1]
////							&& sNeighbor_slave[i]>=segStarts_slave[j])
////				{
////					colIdx[i] = j;
////					assert(edgeNeiSeg_slave[i]==j);
////					colNum[j]++;
////				}
////			}
//		} else
//		{
//			edgeNeiSeg_slave[i]=BLOCKNUM64K;
////			colIdx[i]=BLOCKNUM64K;
////			colNum[BLOCKNUM64K]++;
//		}
//	}
//if(myId==DEBUGID) printArray("%d",edgeNeiSeg_slave,edgeLen);

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
	copyArray(upper_slave,   para_s.frontEdgeData,
				edgeLen, edgeStarts_slave[myId], recvEdges, 0);
	copyArray(diag_slave,    para_s.selfConnData,
				cellLen, segStarts_slave[myId], 0, 1);

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

	volatile swInt *bonus;
	void *tmp;
	int k,dims,idim;
	dims = 0;
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		dims = dims > upper_slave.fArrayDims[i]
			? dims : upper_slave.fArrayDims[i];
	}
	LDM_NEW(tmp,swFloat,edgeLen*dims);
	swFloat *tmpFloat = (swFloat*)tmp;
	swInt   *tmpInt   = (swInt*)tmp;

	LDM_NEW(bonus,int,BLOCKNUM64K+1);
	swInt *faceId;
	LDM_NEW(faceId,swInt,edgeLen);
	for(i=0;i<BLOCKNUM64K+1;i++) {bonus[i]=0;}
	for(i=0;i<edgeLen;i++)
	{
		faceId[accuColNum_slave[edgeNeiSeg_slave[i]+1]] = i;
		accuColNum_slave[edgeNeiSeg_slave[i]+1]++;
	}

	// Reorder Upper
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
//		printf("%d,%f\n",edgeLen,upper_slave.floatArrays[i][0]);
		dims = upper_slave.fArrayDims[i];
		if(edgeLen>0) {
		DMA_Get(tmpFloat,&upper->floatArrays[i][dims*edgeStarts_slave[myId]],edgeLen*dims*sizeof(swFloat));
		}
		for(idim=0;idim<dims;idim++)
		{
			for(j=0;j<edgeLen;j++)
			{
//				if(i==1) printf("%d,%d,%d,%f,%f\n",j,recvEdges,faceId[j],upper_slave.floatArrays[i][(j+recvEdges)*dims+idim],tmpFloat[faceId[j]*dims+idim]);
				upper_slave.floatArrays[i][(j+recvEdges)*dims+idim]
					= tmpFloat[faceId[j]*dims+idim];
			}
		}
	}


	// Reorder neighbor
	if(edgeLen>0) {
    DMA_Get(tmpInt,&neighbor[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
	}
	for(i=0;i<edgeLen;i++)
	{
		neighbor_slave[i+recvEdges] = tmpInt[faceId[i]];
	}

  // Reorder owner
	if(edgeLen>0) {
    DMA_Get(tmpInt,&owner[edgeStarts_slave[myId]],
				edgeLen*sizeof(swInt));
	}
	for(i=0;i<edgeLen;i++)
	{
		owner_slave[i+recvEdges] = tmpInt[faceId[i]];
	}

	// Transform the sorted data to master core
	swInt wrtStart = para_s.wrtStarts[myId+spIndex*BLOCKNUM64K];
	swInt wrtEnd   = para_s.wrtStarts[myId+spIndex*BLOCKNUM64K+1];
	int length = wrtEnd-wrtStart;
	// Output owner/neighbor
	if(length>0)
	{
//		if(myId==7) printf("***%d,%d,%d***\n",segStarts_slave[myId],accuColNum_slave[BLOCKNUM64K],owner_slave[accuColNum_slave[BLOCKNUM64K]]);
		DMA_Put(&rOwner[wrtStart],
				&owner_slave[accuColNum_slave[BLOCKNUM64K]+recvEdges],
				length*sizeof(swInt));
		DMA_Put(&rNeighbor[wrtStart],
				&neighbor_slave[accuColNum_slave[BLOCKNUM64K]+recvEdges],
				length*sizeof(swInt));
	}

////	quickSort(&edgeGlobalIdx[0],edgeLen);
//
////	for(i=0;i<edgeLen;i++) {assert(recvEdges[i]==sNeighbor_slave[i]);}
//
//    swInt ownNeiSendIdx[BLOCKNUM64K],ownNeiRecvIdx[BLOCKNUM64K];
//	swInt ownNeiRecvList[BLOCKNUM64K];
	swInt *ownNeiSendIdx,*ownNeiRecvIdx,*ownNeiRecvList;
	LDM_NEW(ownNeiSendIdx,swInt,BLOCKNUM64K);
	LDM_NEW(ownNeiRecvIdx,swInt,BLOCKNUM64K);
	LDM_NEW(ownNeiRecvList,swInt,BLOCKNUM64K);
	int lastId,startIdx;
//    for(i=0;i<accuColNum_slave[myId+1];i++)
//	{
//		neighbor_slave[i] = sNeighbor_slave[i];
//	}
	swInt fArrayNum     = vertex->fArrayNum; 
	if(fArrayNum>0)
	{
		sVertex_slave.fArrayNum   = fArrayNum;
		sVertex_slave.fArraySizes = cellLen; 
		LDM_NEW(sVertex_slave.fArrayDims, swInt, fArrayNum); 
		LDM_NEW(sVertex_slave.fArrayInOut, swInt, fArrayNum); 
		DMA_Get(sVertex_slave.fArrayDims,vertex->fArrayDims,
					fArrayNum*sizeof(swInt)); 
		DMA_Get(sVertex_slave.fArrayInOut,vertex->fArrayInOut,
					fArrayNum*sizeof(swInt));
		LDM_NEW(sVertex_slave.floatArrays, swFloat*, fArrayNum);
	
		int copyInLen = cellLen+edgeLen-accuColNum_slave[myId+1]+recvEdges;
		for(i=0;i<sVertex_slave.fArrayNum;i++) 
		{ 
			dims = sVertex_slave.fArrayDims[i];
			if(sVertex_slave.fArrayInOut[i]==COPYIN)
			{
				LDM_NEW(sVertex_slave.floatArrays[i],swFloat,
							copyInLen*dims);
				DMA_Get(&sVertex_slave.floatArrays[i][recvEdges*dims], 
					&vertex->floatArrays[i][dims*segStarts_slave[myId]], 
					cellLen*dims*sizeof(swFloat)); 
			} else
			{
				LDM_NEW(sVertex_slave.floatArrays[i],swFloat,cellLen*dims);
				DMA_Get(sVertex_slave.floatArrays[i],
					&vertex->floatArrays[i][dims*segStarts_slave[myId]],
					cellLen*dims*sizeof(swFloat));
			}
		}
	}
	
	// Communicate neighbor
//	if(isLowerExisted)
//	{
		for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
		for(i=0;i<_total_send_pcg;i++)
		{
			if(myId>_sPacks[i].dst_id) continue;
			startIdx
				= accuColNum_slave[_sPacks[i].dst_id]
				+ ownNeiSendIdx[_sPacks[i].dst_id] + recvEdges;
			for(j=0;j<_sPacks[i].cva;j++)
			{
				_sPacks[i].data[j] = neighbor_slave[startIdx+j];
			}
			ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].cva;
		}
		transform_data();
#if 1
		recvNum=0;
		lastId=0;
		k=0;
		for(i=0;i<BLOCKNUM64K;i++) {ownNeiRecvList[i]=0;}
		for(i=0;i<_total_send_pcg;i++)
		{
			if(_rPacks[i].src_id>myId) continue;
			for(j=0;j<_rPacks[i].cva;j++)
			{
				neighbor_slave[recvNum] = _rPacks[i].data[j];
				recvNum++;
			}
			if(_rPacks[i].src_id!=lastId) {k++;lastId=_rPacks[i].src_id;}
			// TODO
			// This will affect communicate X
			ownNeiRecvList[k+1]=recvNum;
			ownNeiRecvIdx[_rPacks[i].src_id]=k;
		}
//		swInt **sdata,**rdata;
//		swInt *ssizes,*rsizes;
//		swInt rnum=0,snum=0;
//		LDM_NEW(sdata,swInt*,BLOCKNUM64K);
//		LDM_NEW(rdata,swInt*,BLOCKNUM64K);
//		LDM_NEW(ssizes,swInt,BLOCKNUM64K);
//		LDM_NEW(rsizes,swInt,BLOCKNUM64K);
//		for(i=0;i<myId;i++)
//		{
//			sdata[i] = NULL;
//			rdata[i] = &neighbor_slave[rnum];
//			ssizes[i] = 0;
//			rsizes[i] = colSegNum_slave[i];
//			rnum += rsizes[i];
//		}
//		sdata[myId]  = NULL;
//		rdata[myId]  = NULL;
//		ssizes[myId] = 0;
//		rsizes[myId] = 0;
//		snum = colSegNum_slave[myId];
//		for(i=myId+1;i<BLOCKNUM64K;i++)
//		{
//			sdata[i] = &neighbor_slave[recvEdges+snum];
//			rdata[i] = NULL;
//			ssizes[i] = colSegNum_slave[i];
//			rsizes[i] = 0;
//			snum += ssizes[i];
//		}
//		groupTrans(sdata,ssizes,rdata,rsizes);
//	}
// ***********************************************************************
// Sparse block computation
// ***********************************************************************
	int iDim;
	upper_slave.fArraySizes
		= accuColNum_slave[BLOCKNUM64K+1]-accuColNum_slave[BLOCKNUM64K];
	diag_slave.fArraySizes   = 0;
	sLower_slave.fArraySizes = 0;
	swFloat* x_M = accessArray(vertex,0);
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		dims = upper_slave.fArrayDims[i];
		upper_slave.floatArrays[i]
			+= (recvEdges+accuColNum_slave[BLOCKNUM64K])*dims;
	}
	Arrays tmpVertex;
	LDM_NEW(tmpVertex.floatArrays, swFloat*, sVertex_slave.fArrayNum);
	LDM_NEW(tmpVertex.fArrayDims,  swFloat,  sVertex_slave.fArrayNum);
	for(i=0;i<sVertex_slave.fArrayNum;i++)
	{
		dims = sVertex_slave.fArrayDims[i];
		if(sVertex_slave.fArrayInOut[i]==COPYIN)
		{
			tmpVertex.floatArrays[i] = vertex->floatArrays[i];
		} else
		{
			tmpVertex.floatArrays[i]
				= sVertex_slave.floatArrays[i]-segStarts_slave[myId]*dims;
		}
	}

	para->fun_slave(&sLower_slave,&upper_slave,&diag_slave,&tmpVertex,
				&owner_slave[accuColNum_slave[BLOCKNUM64K]+recvEdges],
				&neighbor_slave[accuColNum_slave[BLOCKNUM64K]+recvEdges]);

//	swFloat* upperA	= accessArray(&upper_slave, 0);
//	swFloat* bTmp		= accessArray(&tmpVertex, 0);
//	swInt edgeNumber = getArraySize(&upper_slave);
// 	for(iDim=0;iDim<dims;iDim++)
// 	{
//		for(i=0;i<edgeNumber;i++)
// 		{
////if(owner_slave[i+accuColNum_slave[BLOCKNUM64K]+recvEdges]*dims+iDim==43752) printf("sparUpper: %f,%f,%f\n",bTmp[owner_slave[i+accuColNum_slave[BLOCKNUM64K]+recvEdges]*dims+iDim],upperA[i*dims+iDim],x_M[neighbor_slave[i+accuColNum_slave[BLOCKNUM64K]+recvEdges]*dims+iDim]);
// 			bTmp[owner_slave[i+accuColNum_slave[BLOCKNUM64K]+recvEdges]*dims+iDim]
// 				+= upperA[i*dims+iDim];
////				*  x_M[neighbor_slave[i+accuColNum_slave[BLOCKNUM64K]+recvEdges]*dims+iDim];
// 		}
// 	}
//	for(i=0;i<sVertex_slave.fArrayNum;i++)
//	{
//		if(sVertex_slave.fArrayInOut[i]==COPYIN) continue;
//		dims = sVertex_slave.fArrayDims[i];
//		sVertex_slave.floatArrays[i] += segStarts_slave[myId]*dims;
//	}

//	swInt accuConnSeg[2],segStarts_spar[2],edgeStarts_spar[2];
//	accuConnSeg[0] = para_s.accuConnSeg[row];
//	accuConnSeg[1] = para_s.accuConnSeg[row+1];
//	swInt segNum_spar = accuConnSeg[1]-accuConnSeg[0];
//	if(segNum_spar>0)
//	{
//		swInt *endSeg;
//		LDM_NEW(endSeg,swInt,segNum_spar);
//		DMA_Get(endSeg,&para_s.endSeg[accuConnSeg[0]],
//				segNum_spar*sizeof(swInt));
//		for(i=0;i<segNum_spar;i++)
//		{
//if(myId==0)			printf("%d,%d\n",i,endSeg[i]);
////			segStarts_spar[0]  = segStarts[endSeg[i]];
////			segStarts_spar[1]  = segStarts[endSeg[i]+1];
////			edgeStarts_spar[0] = edgeStarts[endSeg[i]];
////			edgeStarts_spar[1] = edgeStarts[endSeg[i]+1];
////			int dims = vertex->fArrayDims[0];
////			int len = segStarts_spar[1]-segStarts_spar[0];
//////			DMA_Get(tmpInt,&vertex->floatArrays[0][segStarts_spar[0]],
//////						len*dims*sizeof(swInt));
//////			for(i=edgeStarts_spar[0];i<edgeStarts_spar[1];i++)
//////			{
//////				neighbor_slave[i]
//		}
//	}
// ***********************************************************************

	// Communicate X
	int iArray;
	swInt idx;
	swInt *ownNeiSendIdxUToL,*ownNeiSendIdxLToU;
	LDM_NEW(ownNeiSendIdxUToL,swInt,BLOCKNUM64K);
	LDM_NEW(ownNeiSendIdxLToU,swInt,BLOCKNUM64K);
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
						+ ownNeiSendIdxUToL[_sPacks[i].dst_id] + recvEdges;
					for(j=0;j<_sPacks[i].cva;j++)
					{
						idx = owner_slave[startIdx+j]
							- segStarts_slave[myId]+recvEdges;
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
						idx = neighbor_slave[startIdx+j]
							- segStarts_slave[myId] + recvEdges;
						_sPacks[i].data[j]=sVertex_slave.floatArrays[iArray][idx*dims+iDim];
//				if(myId==21 && _sPacks[i].dst_id==0) printf("packs: %d,%d,%d,%d,%f,%d\n",idx,recvEdges,neighbor_slave[startIdx+j],startIdx,_sPacks[i].data[j],ownNeiRecvIdx[_sPacks[i].dst_id]);
					}
					ownNeiSendIdxLToU[_sPacks[i].dst_id]+=_sPacks[i].cva;
				}
			}
			transform_data();
			int upVertexIdx  = recvEdges+cellLen;
			int lowVertexIdx = 0;
		    for(i=0;i<_total_send_pcg;i++)
			{
				if(_rPacks[i].src_id>myId)
				{
					for(j=0;j<_rPacks[i].cva;j++)
					{
//if(myId==0 && spIndex==0 && upVertexIdx==396) printf("%d,%d,%d,%f\n",i,j,_rPacks[i].src_id,_rPacks[i].data[j]);
						sVertex_slave.floatArrays[iArray][upVertexIdx*dims+iDim] = _rPacks[i].data[j];
						upVertexIdx++;
					}
				} else
				{
					for(j=0;j<_rPacks[i].cva;j++)
					{
						sVertex_slave.floatArrays[iArray][lowVertexIdx*dims+iDim] = _rPacks[i].data[j];
						lowVertexIdx++;
					}
				}
			}
//		if(upVertexIdx!=copyInLen && myId==DEBUGID) printf("%d,%d,%d,%d\n",upVertexIdx,copyInLen,recvEdges,cellLen);
//		assert(lowVertexIdx==recvEdges);
		}
	}

	for(i=0;i<accuColNum_slave[BLOCKNUM64K]+recvEdges;i++)
	{
		owner_slave[i]     -= segStarts_slave[myId];
		neighbor_slave[i] -= segStarts_slave[myId];
	}
//	for(i=0;i<recvEdges;i++)
//	{
//		neighbor_slave[i]  -= segStarts_slave[myId];
//	}
	for(i=accuColNum_slave[myId+1];i<accuColNum_slave[BLOCKNUM64K];i++)
	{
		neighbor_slave[i+recvEdges] = cellLen+i-accuColNum_slave[myId+1];
	}

// ************************************************************************
// diagonal computation
// ************************************************************************
	swFloat* diagA	= accessArray(&diag_slave, 0);
	swFloat* b		= accessArray(&sVertex_slave, 1);

    sLower_slave.fArraySizes  = 0;
    upper_slave.fArraySizes  = 0;
	diag_slave.fArraySizes   = cellLen;
	for(i=0;i<sVertex_slave.fArrayNum;i++)
	{
		if(sVertex_slave.fArrayInOut[i]==COPYOUT) continue;
		dims = sVertex_slave.fArrayDims[i];
		sVertex_slave.floatArrays[i] += recvEdges*dims;
	}

	swFloat* x		= accessArray(&sVertex_slave, 0);
	para->fun_slave(&sLower_slave,&upper_slave,&diag_slave,
				&sVertex_slave,owner_slave,neighbor_slave);
//	swInt vertexNum = getArraySize(&diag_slave);
//	x = accessArray(&sVertex_slave, 0);
//	dims = diag_slave.fArrayDims[0];
//	swInt ivertex;
//	for( ivertex = 0; ivertex < vertexNum; ivertex++)
//	{
//		for(iDim=0;iDim<dims;iDim++)
//		{
////if((ivertex+segStarts_slave[myId])*dims+iDim==43752) printf("diag:%f,%f,%f\n",b[ivertex*dims+iDim],diagA[ivertex*dims+iDim],x[ivertex*dims+iDim]);
//			b[ivertex*dims+iDim]
//				+=diagA[ivertex*dims+iDim]*x[ivertex*dims+iDim];
//		}
//	}
// ************************************************************************
// diagonal computation
// ************************************************************************


// ************************************************************************
// Compute the upper triangle
// ************************************************************************
	diag_slave.fArraySizes  = 0;
	upper_slave.fArraySizes = accuColNum_slave[BLOCKNUM64K];
	//frontEdge computation
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		dims = upper_slave.fArrayDims[i];
		upper_slave.floatArrays[i] -= accuColNum_slave[BLOCKNUM64K]*dims;
	}
	para->fun_slave(&sLower_slave,&upper_slave,&diag_slave,&sVertex_slave,
				&owner_slave[recvEdges],&neighbor_slave[recvEdges]);

//	edgeNumber = getArraySize(&upper_slave);
//	upperA	= accessArray(&upper_slave, 0);
//	b = accessArray(&sVertex_slave, 1);
//	x = accessArray(&sVertex_slave, 0);
//	swInt iedge;
//	dims = upper_slave.fArrayDims[0];
//	for( iedge = recvEdges; iedge < edgeNumber+recvEdges; iedge++)
//	{
//		for(iDim=0;iDim<dims;iDim++)
//		{
//if((owner_slave[iedge]+segStarts_slave[myId])*dims+iDim==43752) printf("upper:%d,%d,%d,%f,%f,%f\n",iedge,owner_slave[iedge],neighbor_slave[iedge],b[owner_slave[iedge]*dims+iDim],upperA[(iedge-recvEdges)*dims+iDim],x[neighbor_slave[iedge]*dims+iDim]);
//			b[owner_slave[iedge]*dims+iDim]
//				+= upperA[(iedge-recvEdges)*dims+iDim]
//				* x[neighbor_slave[iedge]*dims+iDim];
//		}
//	}
// ************************************************************************
// Compute the upper triangle
// ************************************************************************
	// Transform the sorted data to master core
	// Output upperData
	for(i=0;i<upper_slave.fArrayNum;i++)
	{
		dims = upper_slave.fArrayDims[i];
		upper_slave.floatArrays[i]
			-= recvEdges*dims;
	}
  	for(i=0;i<upper_slave.fArrayNum;i++)
    {
        if(upper_slave.fArrayInOut[i]==COPYIN) continue;
        dims = upper_slave.fArrayDims[i];
        for(idim=0;idim<dims;idim++)
        {
			for(j=0;j<edgeLen;j++)
            {
                tmpFloat[faceId[j]*dims+idim]
                    = upper_slave.floatArrays[i][(j+recvEdges)*dims+idim];
            }
        }
		if(edgeLen>0) {
        DMA_Put(&upper->floatArrays[i][dims*edgeStarts_slave[myId]],
            tmpFloat,edgeLen*dims*sizeof(swFloat));
		}
		destroyRegisterPacks(rlmpi_info);
		return;
	}

//	for(i=0;i<upper_slave.fArrayNum;i++)
//	{
//		if(upper_slave.fArrayInOut[i]==COPYIN) continue;
//		DMA_Put(upper->floatArrays[i],
//			&upper_slave.floatArrays[i]
//			[recvEdges*upper_slave.fArrayDims[i]],
//			edgeLen*upper_slave.fArrayDims[i]*sizeof(swFloat));
//	}

// ************************************************************************
// Lower part
// ************************************************************************
//	copyArray(sLower_slave,  lower,
//				edgeLen, edgeStarts_slave[myId], recvEdges, 0);

	sLower_slave = upper_slave;
	// Reorder lower
	sLower_slave.fArrayNum = lower->fArrayNum;
	sLower_slave.fArraySizes = lower->fArraySizes;
	if(sLower_slave.fArrayNum>0)
	{
		for(i=0;i<sLower_slave.fArrayNum;i++)
		{
			dims = sLower_slave.fArrayDims[i];
			if(edgeLen*dims>0) {
			DMA_Get(tmpFloat,&lower->floatArrays[i][dims*edgeStarts_slave[myId]],edgeLen*dims*sizeof(swFloat));
			}
	//		dims=1;
	//		idim=0;
			for(idim=0;idim<dims;idim++)
			{
				for(j=0;j<edgeLen;j++)
				{
					sLower_slave.floatArrays[i][(j+recvEdges)*dims+idim]
						= tmpFloat[faceId[j]*dims+idim];
				}
			}
		}
	}

	// Communicate Lower
	if(isLowerExisted)
	{
	    for(iArray=0;iArray<sLower_slave.fArrayNum;iArray++)
		{
			dims = sLower_slave.fArrayDims[iArray];
			for(iDim=0;iDim<dims;iDim++)
			{
				for(i=0;i<BLOCKNUM64K;i++) {ownNeiSendIdx[i]=0;}
				for(i=0;i<_total_send_pcg;i++)
				{
					if(myId>_sPacks[i].dst_id) continue;
					startIdx = accuColNum_slave[_sPacks[i].dst_id]
						+ ownNeiSendIdx[_sPacks[i].dst_id] + recvEdges;
					for(j=0;j<_sPacks[i].cva;j++)
					{
						_sPacks[i].data[j] = sLower_slave.floatArrays[iArray][(startIdx+j)*dims+iDim];
					}
					ownNeiSendIdx[_sPacks[i].dst_id]+=_sPacks[i].cva;
				}
				transform_data();
				startIdx=0;
				for(i=0;i<_total_send_pcg;i++)
				{
					if(_rPacks[i].src_id>myId) continue;
					for(j=0;j<_rPacks[i].cva;j++)
					{
						sLower_slave.floatArrays[iArray][startIdx*dims+iDim]
							= _rPacks[i].data[j];
						startIdx++;
					}
				}
			}
		}
	}
//if(myId==0) printArray("%f",upVertex_slave.floatArrays[0],cellLen+edgeLen-colNum[myId]);
//if(myId==1) printArray("%d",sNeighbor_slave,edgeLen);
//if(myId==1) printArray("%f",upVertex_slave.floatArrays[0],edgeLen);


//if(myId==4) printArray("%f",upVertex_slave.floatArrays[1],cellLen);

	// Compute the lower triangle
	for(i=0;i<recvEdges;i++)
	{
		owner_slave[i] = i;
	}
	for(i=0;i<accuColNum_slave[myId+1];i++)
	{
		owner_slave[i+recvEdges] += recvEdges;
	}

// ************************************************************************
// backEdge computation
// ************************************************************************
    upper_slave.fArraySizes  = 0;
	diag_slave.fArraySizes   = 0;
	sLower_slave.fArraySizes = recvEdges+accuColNum_slave[myId+1];
	for(i=0;i<sVertex_slave.fArrayNum;i++)
	{
		if(sVertex_slave.fArrayInOut[i]==COPYOUT) continue;
		dims = sVertex_slave.fArrayDims[i];
		sVertex_slave.floatArrays[i] -= recvEdges*dims;
	}

//	x		= accessArray(&sVertex_slave, 0);
//	b		= accessArray(&sVertex_slave, 1);
	para_s.fun_slave(&sLower_slave,&upper_slave,&diag_slave,
				&sVertex_slave,owner_slave,neighbor_slave);
	//backEdge computation
//	swFloat* lowerA	= accessArray(&sLower_slave, 0);
//	dims = sLower_slave.fArrayDims[0];
//	edgeNumber = getArraySize(&sLower_slave);
//	for(iDim=0;iDim<dims;iDim++)
//	{
//		for( iedge = 0; iedge < edgeNumber; iedge++)
//		{
//if((neighbor_slave[iedge]+segStarts_slave[myId])*dims+iDim==43752) printf("lower:%d,%d,%d,%d,%f,%f,%f\n",iedge,iedge*dims+iDim,owner_slave[iedge],neighbor_slave[iedge],b[neighbor_slave[iedge]*dims+iDim],lowerA[iedge*dims+iDim],x[owner_slave[iedge]*dims+iDim]);
//			b[neighbor_slave[iedge]*dims+iDim]
//				+= lowerA[iedge*dims+iDim]*x[owner_slave[iedge]*dims+iDim];
//		}
//	}
	// Sparse block computation
//	swFloat* b_M = accessArray(vertex,1);
//	swFloat tmpVal;
//	for(i=accuColNum_slave[BLOCKNUM64K]+recvEdges;
//				i<accuColNum_slave[BLOCKNUM64K+1]+recvEdges;i++)
//	{
//		for(iDim=0;iDim<dims;iDim++)
//		{
////if((neighbor_slave[i])*dims+iDim==26604 && myId==63) printf("lowerSparse:%d,%d,%d,%f,%f,%f,%d\n",i*dims+iDim,owner_slave[i],neighbor_slave[i],b_M[neighbor_slave[i]*dims+iDim],lowerA[i*dims+iDim],x[(owner_slave[i]-segStarts_slave[myId]+recvEdges)*dims+iDim],para_s.segEdgeNum[0]);
//			b_M[neighbor_slave[i]*dims+iDim] +=
//				lowerA[i*dims+iDim] * x[(owner_slave[i]-segStarts_slave[myId]+recvEdges)*dims+iDim];
////if((neighbor_slave[i])*dims+iDim==26604 && myId==63) printf("lowerSparse:%d,%d,%d,%f,%f,%f,%d\n",i*dims+iDim,owner_slave[i],neighbor_slave[i],b_M[neighbor_slave[i]*dims+iDim],lowerA[i*dims+iDim],x[(owner_slave[i]-segStarts_slave[myId]+recvEdges)*dims+iDim],para_s.segEdgeNum[0]);
//		}
//	}

// ************************************************************************
// backEdge computation
// ************************************************************************

	// Combine lower vertex and upper vertex
//	for(i=0;i<lowVertex_slave.fArrayNum;i++)
//	{
//		dims = lowVertex_slave.fArrayDims[i];
//		if(lowVertex_slave.fArrayInOut[i]==COPYOUT)
//		{
//			for(j=0;j<cellLen*dims;j++)
//			{
//				lowVertex_slave.floatArrays[i][j] 
//					+= upVertex_slave.floatArrays[i][j];
//			}
//		}
//	}

//if(myId==0) printArray("%f",lowVertex_slave.floatArrays[1],cellLen);
	// TODO
	// Assume that lower data won't be output
//if(myId==7) printf("%d,%d\n",segStart,segStarts_slave[myId]);
//	lowVertex_slave.fArraySizes   = cellLen;
	diag_slave.fArraySizes        = cellLen;
	upper_slave.fArraySizes       = edgeLen;
	// Output vertexData
	for(i=0;i<sVertex_slave.fArrayNum;i++)
	{
		if(sVertex_slave.fArrayInOut[i]==COPYIN) continue;
		DMA_Put(&vertex->floatArrays[i]
				[segStarts_slave[myId]*sVertex_slave.fArrayDims[i]],
				sVertex_slave.floatArrays[i],
				cellLen*sVertex_slave.fArrayDims[i]*sizeof(swFloat));
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
	// Output lowerData
	if(length>0)
	{
		for(i=0;i<sLower_slave.fArrayNum;i++)
		{
			DMA_Put(&rLower->floatArrays[i]
					[wrtStart*sLower_slave.fArrayDims[i]],
					&sLower_slave.floatArrays[i]
					[(accuColNum_slave[BLOCKNUM64K]+recvEdges)
					*sLower_slave.fArrayDims[i]],
					length*sLower_slave.fArrayDims[i]*sizeof(swFloat));
		}
	}

#endif
	destroyRegisterPacks(rlmpi_info);
//	updt_addl(1, &para_s.segEdgeNum[0]);
//	para_s.segEdgeNum[0]+=1;
//	AO_SET(&para_s.segEdgeNum[0],1);
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


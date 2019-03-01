#include "swMacro.h"
#include "slave.h"
#include <assert.h>
#include "rlmpi.h"
#include "slaveUtils.h"
#include "multiLevelBlockIterator.h"


#define copyArrayWithoutData(srcArray, dstArray, size) \
{ \
	int iArray; \
	(dstArray).fArraySizes = (size); \
	(dstArray).fArrayNum   = (srcArray)->fArrayNum; \
	LDM_NEW((dstArray).fArrayInOut,swInt,(dstArray).fArrayNum); \
	LDM_NEW((dstArray).fArrayDims, swInt,(dstArray).fArrayNum); \
	LDM_NEW((dstArray).floatArrays,swFloat*,(dstArray).fArrayNum); \
	for(iArray=0;iArray<(dstArray).fArrayNum;iArray++) \
	{ \
		(dstArray).fArrayInOut[iArray] \
			= (srcArray)->fArrayInOut[iArray]; \
		(dstArray).fArrayDims[iArray] \
			= (srcArray)->fArrayDims[iArray]; \
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

#define RLC_vertexData(vertexData_s) \
{ \
	swInt *sendIdxRev; \
	swInt idx; \
	LDM_NEW(sendIdxRev, swInt, BLOCKNUM64K); \
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++) \
	{ \
		if(vertexData_s.fArrayInOut[iArray]==COPYOUT) continue; \
		dims = vertexData_s.fArrayDims[iArray]; \
		swFloat *x=vertexData_s.floatArrays[iArray]; \
		for(iDim=0;iDim<dims;iDim++) \
		{ \
			for(j=0;j<BLOCKNUM64K;j++) {sendIdx[j]=0; sendIdxRev[j]=0;} \
			for(ipcg=0;ipcg<_total_send_pcg;ipcg++) \
			{ \
				if(_MYID>_sPacks[ipcg].dst_id) \
				{ \
					edgeIdx = recvList[recvIdx[_sPacks[ipcg].dst_id]] \
						+ sendIdxRev[_sPacks[ipcg].dst_id]; \
					for(j=0;j<_sPacks[ipcg].cva;j++) \
					{ \
						idx = neighbor[edgeIdx+j]-cellStarts_s[_MYID] \
							+ recvEdges; \
						_sPacks[ipcg].data[j] = x[idx*dims+iDim]; \
					} \
					sendIdxRev[_sPacks[ipcg].dst_id]+=_sPacks[ipcg].cva; \
				} else if(_MYID<_sPacks[ipcg].dst_id) \
				{ \
					edgeIdx \
						= blockStarts_s[4*(_sPacks[ipcg].dst_id-_MYID)+2] \
						+ sendIdx[_sPacks[ipcg].dst_id] + recvEdges \
						- blockStarts_s[2]; \
					for(j=0;j<_sPacks[ipcg].cva;j++) \
					{ \
						idx = owner[edgeIdx+j]-cellStarts_s[_MYID] \
							+ recvEdges; \
						_sPacks[ipcg].data[j] = x[idx*dims+iDim]; \
					} \
					sendIdx[_sPacks[ipcg].dst_id]+=_sPacks[ipcg].cva; \
				} \
			} \
			transform_data(); \
			recvNum = recvEdges+cellLen; \
			int recvNumRev = 0; \
			for(ipcg=0;ipcg<_total_recv_pcg;ipcg++) \
			{ \
				if(_rPacks[ipcg].src_id>_MYID) \
				{ \
					for(j=0;j<_rPacks[ipcg].cva;j++) \
					{ \
						x[recvNum*dims+iDim] = _rPacks[ipcg].data[j]; \
						recvNum++; \
					} \
				} else \
				{ \
					for(j=0;j<_rPacks[ipcg].cva;j++) \
					{ \
						x[recvNumRev*dims+iDim] = _rPacks[ipcg].data[j]; \
						recvNumRev++; \
					} \
				} \
			} \
		} \
	} \
} 

#define RLC_backEdgeData(backEdgeData_s) \
{ \
	for(iArray=0;iArray<backEdgeData_s.fArrayNum;iArray++) \
	{ \
		dims = backEdgeData_s.fArrayDims[iArray]; \
		swFloat* lower = backEdgeData_s.floatArrays[iArray]; \
		for(iDim=0;iDim<dims;iDim++) \
		{ \
			for(j=0;j<BLOCKNUM64K;j++) { sendIdx[j]=0; } \
			for(ipcg=0;ipcg<_total_send_pcg;ipcg++) \
			{ \
				if(_MYID>=_sPacks[ipcg].dst_id) continue; \
				edgeIdx = blockStarts_s[4*(_sPacks[ipcg].dst_id-_MYID)+2] \
					+sendIdx[_sPacks[ipcg].dst_id]+recvEdges-blockStarts_s[2]; \
				for(j=0;j<_sPacks[ipcg].cva;j++) \
				{ \
					_sPacks[ipcg].data[j] = lower[(edgeIdx+j)*dims+iDim]; \
				} \
				sendIdx[_sPacks[ipcg].dst_id] += _sPacks[ipcg].cva; \
			} \
			transform_data(); \
			recvNum=0; \
			for(ipcg=0;ipcg<_total_recv_pcg;ipcg++) \
			{ \
				if(_MYID<=_rPacks[ipcg].src_id) continue; \
				for(j=0;j<_rPacks[ipcg].cva;j++) \
				{ \
					lower[recvNum*dims+iDim] = _rPacks[ipcg].data[j]; \
					recvNum++; \
				} \
			} \
		} \
	} \
} 

void multiLevelBlockIterator_e2v_slave(MLB_edge2VertexPara *para)
{
	MLB_edge2VertexPara para_s;
	DMA_Get(&para_s, para, sizeof(MLB_edge2VertexPara));

	INIT_LDM_SPACE(48000);
	
	swInt *cellStarts_s, *blockStarts_s,
		  *owner_s, *neighbor_s;
	LDM_NEW(cellStarts_s,  swInt, BLOCKNUM64K+1);
	LDM_NEW(blockStarts_s, swInt, BLOCKNUM64K*4);

	swInt spIndex = para_s.spIndex;
	swInt mshBlockNum = para_s.mshBlockNum;
	swInt maxEdges = para_s.maxEdges;
	swInt maxCells = para_s.maxCells;
	swInt recvEdges = para_s.recvStarts[_MYID+spIndex*BLOCKNUM64K];
//	recvEdges = 0;

	swInt rowId = _MYID+spIndex*BLOCKNUM64K;
	swInt cpeBlockNum = mshBlockNum*BLOCKNUM64K;
	swInt startBlockIdx = rowId*(1+2*cpeBlockNum-rowId)/2;

	DMA_Get(blockStarts_s, &para_s.blockStarts[4*startBlockIdx],
				BLOCKNUM64K*4*sizeof(swInt));
	DMA_Get(cellStarts_s, &para_s.cellStarts[spIndex*BLOCKNUM64K],
				(BLOCKNUM64K+1)*sizeof(swInt));

	RlmpiInfo *rlmpi_info = &para_s.schedule_data[spIndex];
	initRegisterPacks(rlmpi_info);

	Arrays backEdgeData_s,frontEdgeData_s,vertexData_s,selfConnData_s;

	int iArray;
	swInt *owner,*neighbor;
	int startIdx = blockStarts_s[2];
	int edgeLen = blockStarts_s[4*(BLOCKNUM64K-_MYID-1)+3]-startIdx;
	int localEdgeLen = blockStarts_s[3]-blockStarts_s[2];
	int cellLen = cellStarts_s[_MYID+1]-cellStarts_s[_MYID];
	LDM_NEW(owner, swInt, edgeLen+recvEdges);
	LDM_NEW(neighbor, swInt, edgeLen+recvEdges);

	// Construct array from master core without data copying
	copyArrayWithoutData(para_s.vertexData,    vertexData_s,   cellLen);
	copyArrayWithoutData(para_s.selfConnData,  selfConnData_s, cellLen);
	copyArrayWithoutData(para_s.backEdgeData,  backEdgeData_s,
				edgeLen+recvEdges);
	copyArrayWithoutData(para_s.frontEdgeData, frontEdgeData_s,
				edgeLen+recvEdges);

	backEdgeData_s.fArraySizes = edgeLen;
	frontEdgeData_s.fArraySizes = edgeLen;
	if(edgeLen>0)
	{
		DMA_Get(&owner[recvEdges], &para_s.owner[startIdx],
					edgeLen*sizeof(swInt));
		DMA_Get(&neighbor[recvEdges], &para_s.neighbor[startIdx],
					edgeLen*sizeof(swInt));
		for(iArray=0;iArray<frontEdgeData_s.fArrayNum;iArray++)
		{
			int dims=frontEdgeData_s.fArrayDims[iArray]; 
			LDM_NEW(frontEdgeData_s.floatArrays[iArray],swFloat,
					(recvEdges+edgeLen)*dims);
			DMA_Get(&frontEdgeData_s.floatArrays[iArray][recvEdges*dims],
					&para_s.frontEdgeData->floatArrays[iArray][startIdx*dims],
					edgeLen*dims*sizeof(swFloat));
		}
	}

	// Get vertex data through DMA
	int iedge,dims,iDim;
	int RL_length = cellLen+recvEdges+edgeLen-localEdgeLen;
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		if(vertexData_s.fArrayInOut[iArray]==COPYOUT)
		{
			LDM_NEW(vertexData_s.floatArrays[iArray],swFloat,cellLen*dims);
			DMA_Get(vertexData_s.floatArrays[iArray],
					&para_s.vertexData->floatArrays
					[iArray][cellStarts_s[_MYID]*dims],
					cellLen*dims*sizeof(swFloat));
		} else if(vertexData_s.fArrayInOut[iArray]==COPYIN)
		{
			LDM_NEW(vertexData_s.floatArrays[iArray],swFloat,
					RL_length*dims);
			DMA_Get(&vertexData_s.floatArrays[iArray][recvEdges*dims],
					&para_s.vertexData->floatArrays
					[iArray][cellStarts_s[_MYID]*dims],
					cellLen*dims*sizeof(swFloat));
		}
	}
	
	// Get self connection data through DMA
	for(iArray=0;iArray<selfConnData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&selfConnData_s,iArray);
		LDM_NEW(selfConnData_s.floatArrays[iArray],swFloat,cellLen*dims);
		DMA_Get(selfConnData_s.floatArrays[iArray],
					&para_s.selfConnData->floatArrays
					[iArray][cellStarts_s[_MYID]*dims],
					cellLen*dims*sizeof(swFloat));
	}


	// Register communication: neighbor
	int ipcg,j,recvNum,edgeIdx;
	swInt *sendIdx,*recvIdx,*recvList;
	LDM_NEW(sendIdx,  swInt, BLOCKNUM64K);
	LDM_NEW(recvIdx,  swInt, BLOCKNUM64K);
	LDM_NEW(recvList, swInt, BLOCKNUM64K);
	for(j=0;j<BLOCKNUM64K;j++) { sendIdx[j]=0; }
	for(ipcg=0;ipcg<_total_send_pcg;ipcg++)
	{
		if(_MYID>=_sPacks[ipcg].dst_id) continue;
		edgeIdx = blockStarts_s[4*(_sPacks[ipcg].dst_id-_MYID)+2]
			+ sendIdx[_sPacks[ipcg].dst_id]+recvEdges - blockStarts_s[2];
		for(j=0;j<_sPacks[ipcg].cva;j++)
		{
			_sPacks[ipcg].data[j] = neighbor[edgeIdx+j];
//if(neighbor[startIdx+j]==125793) printf("%d,%d,%f\n",_MYID,startIdx+j-recvEdges+blockStarts_s[2],frontEdgeData_s.floatArrays[0][startIdx+j]);
		}
		sendIdx[_sPacks[ipcg].dst_id] += _sPacks[ipcg].cva;
	}
	transform_data();
	recvNum=0;
	int k=0,lastId=0;
	for(j=0;j<BLOCKNUM64K;j++) { recvList[j]=0; }
	for(ipcg=0;ipcg<_total_recv_pcg;ipcg++)
	{
		if(_MYID<=_rPacks[ipcg].src_id) continue;
		for(j=0;j<_rPacks[ipcg].cva;j++)
		{
			neighbor[recvNum++] = _rPacks[ipcg].data[j];
//if(neighbor[recvNum-1]==128330) printf("%d,%d,%f\n",_MYID,recvNum-1+blockStarts_s[2],backEdgeData_s.floatArrays[0][recvNum-1]);
		}
		if(_rPacks[ipcg].src_id!=lastId){k++;lastId=_rPacks[ipcg].src_id;}
		recvList[k+1]=recvNum;
		recvIdx[_rPacks[ipcg].src_id]=k;
//		recvNum += _rPacks[ipcg].cva;
	}

	// Register communication: vertex data(COPYIN)
	RLC_vertexData(vertexData_s);

	// localization
	for(iedge=recvEdges+localEdgeLen;iedge<recvEdges+edgeLen;iedge++)
	{
		neighbor[iedge] = iedge-localEdgeLen+cellLen+cellStarts_s[_MYID]-recvEdges;
	}

	// Compute Diag
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		if(vertexData_s.fArrayInOut[iArray]==COPYIN)
		{
			vertexData_s.floatArrays[iArray] += recvEdges*dims;
		}
	}
	frontEdgeData_s.fArraySizes = 0;
	backEdgeData_s.fArraySizes = 0;
	selfConnData_s.fArraySizes = cellLen;
	para_s.fun_slave(&backEdgeData_s, &frontEdgeData_s, &selfConnData_s,
				&vertexData_s, &owner[recvEdges], &neighbor[recvEdges]);

	// Compute Upper
	Arrays tmpVertexData;
	LDM_NEW(tmpVertexData.floatArrays, swFloat*, vertexData_s.fArrayNum);
	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		tmpVertexData.floatArrays[iArray] 
			= vertexData_s.floatArrays[iArray]-cellStarts_s[_MYID]*dims;
	}
	for(iArray=0;iArray<frontEdgeData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&frontEdgeData_s,iArray);
		frontEdgeData_s.floatArrays[iArray] += recvEdges*dims;
	}
	frontEdgeData_s.fArraySizes = edgeLen;
	backEdgeData_s.fArraySizes = 0;
	selfConnData_s.fArraySizes = 0;
	para_s.fun_slave(&backEdgeData_s, &frontEdgeData_s, &selfConnData_s,
				&tmpVertexData, &owner[recvEdges], &neighbor[recvEdges]);


	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		if(vertexData_s.fArrayInOut[iArray]==COPYIN)
		{
			tmpVertexData.floatArrays[iArray] -= recvEdges*dims;
		}
	}

#if 1
//	swFloat *b = accessArray(&tmpVertexData,0);
//	swFloat *x = accessArray(&tmpVertexData,1);
//	swFloat *upper = accessArray(&frontEdgeData_s,0);
//	dims = getArrayDims(&frontEdgeData_s,0);
////	b -= cellStarts_s[_MYID];
//	for(iedge=0;iedge<edgeLen;iedge++)				
//	{
//		for(iDim=0;iDim<dims;iDim++)
//		{
//if(owner[iedge+recvEdges]==114912) printf("upper: %d,%d,%d,%f,%f,%f\n",_MYID,neighbor[iedge+recvEdges],recvEdges+cellStarts_s[_MYID],b[owner[iedge+recvEdges]*dims+iDim],upper[iedge*dims+iDim],x[neighbor[iedge+recvEdges]*dims+iDim]);
////			printf("%d,%d,%d,%d,%d,%d,%d\n",iedge,edgeLen,neighbor[iedge],cellStarts_s[iblk+_MYID+1],maxCells,ownCellLen,neiCellLen);
//			b[(owner[iedge+recvEdges])*dims+iDim]
//				+= upper[(iedge)*dims+iDim]*x[neighbor[iedge+recvEdges]*dims+iDim];
//		}
//	}


// ***********************************************************************
// Backward edge
// ***********************************************************************
	if(edgeLen>0)
	{
		for(iArray=0;iArray<backEdgeData_s.fArrayNum;iArray++)
		{
			int dims=backEdgeData_s.fArrayDims[iArray]; 
			backEdgeData_s.floatArrays[iArray]
				= frontEdgeData_s.floatArrays[iArray]-recvEdges*dims;
			DMA_Get(&backEdgeData_s.floatArrays[iArray][recvEdges*dims],
					&para_s.backEdgeData->floatArrays[iArray][startIdx*dims],
					edgeLen*dims*sizeof(swFloat));
		}
	}

	// Register communication: backEdgeData
	RLC_backEdgeData(backEdgeData_s);

	// localization
	for(iedge=recvEdges;iedge<recvEdges+localEdgeLen;iedge++)
	{
		owner[iedge] += recvEdges;
	}
	for(iedge=0;iedge<recvEdges;iedge++)
	{
		owner[iedge] = iedge+cellStarts_s[_MYID];
	}

	// Compute Lower
	backEdgeData_s.fArraySizes = localEdgeLen+recvEdges;
	selfConnData_s.fArraySizes = 0;
	frontEdgeData_s.fArraySizes = 0;
	para_s.fun_slave(&backEdgeData_s, &frontEdgeData_s, &selfConnData_s,
				&tmpVertexData, owner, neighbor);

//	b = accessArray(&tmpVertexData,0);
//	x = accessArray(&tmpVertexData,1);
//	swFloat *lower = accessArray(&backEdgeData_s,0);
//	dims = getArrayDims(&backEdgeData_s,0);
////	b -= cellStarts_s[_MYID];
//	swInt length = localEdgeLen+recvEdges;
//	for(iedge=0;iedge<length;iedge++)				
//	{
//		for(iDim=0;iDim<dims;iDim++)
//		{
//if(neighbor[iedge*dims+iDim]==114912) printf("lower: %d,%d,%f,%f,%f\n",iedge,blockStarts_s[2],b[neighbor[iedge*dims+iDim]],lower[iedge*dims+iDim],x[owner[iedge]*dims+iDim]);
////			printf("%d,%d,%d,%d,%d,%d,%d\n",iedge,edgeLen,neighbor[iedge],cellStarts_s[iblk+_MYID+1],maxCells,ownCellLen,neiCellLen);
//			b[(neighbor[iedge])*dims+iDim]
//				+= lower[iedge*dims+iDim]*x[owner[iedge]*dims+iDim];
//		}
//	}

	for(iArray=0;iArray<vertexData_s.fArrayNum;iArray++)
	{
		dims = getArrayDims(&vertexData_s,iArray);
		if(vertexData_s.fArrayInOut[iArray]==COPYOUT)
		{
			DMA_Put(&para_s.vertexData->floatArrays
					[iArray][cellStarts_s[_MYID]*dims],
					vertexData_s.floatArrays[iArray],
					cellLen*dims*sizeof(swFloat));
		}
	}

#endif
	destroyRegisterPacks(rlmpi_info);
}

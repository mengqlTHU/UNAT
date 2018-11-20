#include <stdlib.h>
#include <stdio.h>
#include "slave.h"
#include "simd.h"
#include "dma.h"
#include "rlmpi.h"
#include "swMacro.h"
#include "iterator_struct.h"

__thread_local volatile int myId;

extern void (*operatorFunPointer_host)(Arrays*, Arrays*, swInt*, swInt*);
extern swInt spIndex, maxXNum, maxCell, maxEdge, cpeBlockNum;
extern swInt *blockStarts, *cellStarts, *ownNeiSendList, *owner, *neighbor;
extern swFloat *upper, *lower, *diag, *x, *b;

void func()
{
	myId = athread_get_id(-1);
	if(myId==0) printf("edge2VertexIteration_slave\n");
	swInt i,j,k,row,col,rowId;

	volatile swInt total_send_pcg = _table_ldm.nPUTR + _table_ldm.nPutrSameRow + _table_ldm.nPutcSameCol;
	volatile swInt totalLength
		= maxXNum*3*sizeof(swFloat)
		+ maxEdge*2*(sizeof(swFloat)+sizeof(swInt))
		+ maxCell*2*sizeof(swInt)
		+ maxCell*sizeof(swFloat)
		+ BLOCKNUM64K*8*sizeof(swInt)
		+ total_send_pcg*2*32;
	if(totalLength > 64*1024*0.95 && myId==0){
//	if(myId == 0){
		printf("The allocated LDM exceeds 64KB, the memory size is %d\n",
					totalLength);
	}
	else if(totalLength < 64*1024*0.8 && myId==0){
		printf("The allocated LDM is too small, the memory size is %d\n",
					totalLength);
	}

	volatile swFloat sendX_slave[maxXNum],recvX_slave[maxXNum],dmaX_slave;
	volatile swFloat b_slave[maxXNum];
	volatile swFloat diagUpper_slave[maxCell];
	volatile swInt diagOwner_slave[maxCell],diagNeighbor_slave[maxCell];

	volatile swInt cell_slave[BLOCKNUM64K+1];

	volatile swInt blockStarts_slave[BLOCKNUM64K*4];
    volatile swInt startIdx,blockLen,startBlockIdx,endIdx;
	volatile swInt swIntegral,remainder,kk,l,kkk,jj;
	
	//寄存器通信发送面数据
	volatile swFloat sLower_slave[maxEdge],rLower_slave[maxEdge];
	volatile swInt sNeighbor_slave[maxEdge];
	volatile swInt rNeighbor_slave[maxEdge];

	volatile swInt sPackIdx,testNum,nonZeroNum;
	volatile swInt ownNeiRecvList_slave[BLOCKNUM64K];
	volatile swInt ownNeiRecvIdx[BLOCKNUM64K];
	volatile swInt ownNeiSendIdx[BLOCKNUM64K];


//	sPackIdx = 0;
//	row = myId;
//	row = myId + index*BLOCKNUM64K;
//	startBlockIdx = row*(1+2*cpeBlockNums-row)/2;
//	cell_slave[0] = cellStarts[row];
//	cell_slave[1] = cellStarts[row+1];
////	DMA_Get(&cell_slave[0],&cellStarts[row],2*sizeof(swInt));
//	DMA_Get(&blockStarts_slave[0],&blockStarts[4*startBlockIdx],
//				BLOCKNUM64K*4*sizeof(swInt));
//	DMA_Get(&b_slave[0],&b_MLBSlave[cell_slave[0]],
//				(cell_slave[1]-cell_slave[0])*sizeof(swFloat));
//
//
//	startIdx = blockStarts_slave[6];
//	blockLen = myId==BLOCKNUM64K-1 ?
//		0 : blockStarts_slave[4*(BLOCKNUM64K-myId-1)+3]-startIdx;
//	if(blockLen>0){
//		DMA_Get(&sLower_slave[0],&lowerReordered[startIdx],
//					blockLen*sizeof(swFloat));
//		DMA_Get(&sNeighbor_slave[0],&ownReordered[startIdx],
//					blockLen*sizeof(swFloat));
//	}
//
//	DMA_Get(&sendX_slave[0],&xReordered[cell_slave[0]],
//				(cell_slave[1]-cell_slave[0])*sizeof(swFloat));
//
//
//	//communicate Lower
//	for(j=0;j<BLOCKNUM64K;j++){ownNeiSendIdx[j]=0;}
//	for(j=0;j<total_send_pcg;j++){
//		if(myId>_sPacks[j].dst_id) continue;
//		startIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+2]
//			+ ownNeiSendIdx[_sPacks[j].dst_id];
//		endIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+3];
//		kk=startIdx-blockStarts_slave[6];
//		if(endIdx-startIdx>=dataLen){
//			for(i=0;i<dataLen;i++){
//				_sPacks[j].data[i]
//					= sLower_slave[kk+i]
//					* sendX_slave[sNeighbor_slave[kk+i]-cell_slave[0]];
//			}
//			ownNeiSendIdx[_sPacks[j].dst_id]+=dataLen;
//			_sPacks[j].indM = dataLen;
//		}else{
//			for(i=0;i<endIdx-startIdx;i++){
//				_sPacks[j].data[i]
//					= sLower_slave[kk+i]
//					* sendX_slave[sNeighbor_slave[kk+i]-cell_slave[0]];
//			}
//			ownNeiSendIdx[_sPacks[j].dst_id]+=endIdx-startIdx;
//			_sPacks[j].indM = endIdx-startIdx;
//		}
//	}
//
//	transform_data();
//	nonZeroNum=0;
//	for(i=0;i<total_send_pcg;i++){
//		if(_rPacks[i].src_id<myId){
//			for(j=0;j<_rPacks[i].indM;j++){
//				rLower_slave[nonZeroNum] = _rPacks[i].data[j];
//				nonZeroNum++;
//			}
//		}
//	}
//	if(blockLen>0){
//		startIdx = blockStarts_slave[6];
//		DMA_Get(&sNeighbor_slave[0],&neiReordered[startIdx],
//					blockLen*sizeof(swInt));
//	}
//
//	//communicate Neighbor
//	for(j=0;j<BLOCKNUM64K;j++){ownNeiSendIdx[j]=0;}
//	for(j=0;j<total_send_pcg;j++){
//		if(myId>_sPacks[j].dst_id) continue;
//		startIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+2]
//			+ ownNeiSendIdx[_sPacks[j].dst_id];
//		endIdx = blockStarts_slave[4*(_sPacks[j].dst_id-myId)+3];
//		kk=startIdx-blockStarts_slave[6];
//		if(endIdx-startIdx>=dataLen){
//			for(i=0;i<dataLen;i++){
//				_sPacks[j].data[i]=sNeighbor_slave[kk+i];
//			}
//			ownNeiSendIdx[_sPacks[j].dst_id]+=dataLen;
//			_sPacks[j].indM = dataLen;
//		}else{
//			for(i=0;i<endIdx-startIdx;i++){
//				_sPacks[j].data[i]=sNeighbor_slave[kk+i];
//			}
//			ownNeiSendIdx[_sPacks[j].dst_id]+=endIdx-startIdx;
//			_sPacks[j].indM = endIdx-startIdx;
//		}
//	}
//
//
//	transform_data();
//	nonZeroNum=0;
//	swInt lastId=0;
//	k=0;
//	for(i=0;i<BLOCKNUM64K;i++){ownNeiRecvList_slave[i] = 0;}
//	for(i=0;i<total_send_pcg;i++){
//		if(_rPacks[i].src_id<myId){
//			for(j=0;j<_rPacks[i].indM;j++){
//				rNeighbor_slave[nonZeroNum] = _rPacks[i].data[j];
//				nonZeroNum++;
//			}
//			if(_rPacks[i].src_id!=lastId){k++;lastId=_rPacks[i].src_id;}
//			ownNeiRecvList_slave[k+1]=nonZeroNum;
//			ownNeiRecvIdx[_rPacks[i].src_id]=k;
//		}
//	}
//
//	//计算稠密块中下三角
//	for(i=0;i<nonZeroNum;i++){
//		b_slave[rNeighbor_slave[i]-cell_slave[0]] += rLower_slave[i];
//	}
//
//	//communicate X
//	for(j=0;j<BLOCKNUM64K;j++){ownNeiSendIdx[j]=0;}
//	for(j=0;j<total_send_pcg;j++){
//		if(myId<_sPacks[j].dst_id) continue;
//		startIdx = ownNeiRecvList_slave[ownNeiRecvIdx[_sPacks[j].dst_id]]
//			+ ownNeiSendIdx[_sPacks[j].dst_id];
//		endIdx = ownNeiRecvList_slave[ownNeiRecvIdx[_sPacks[j].dst_id]+1];
//		kk=startIdx;
//		if(endIdx-startIdx>=dataLen){
//			for(i=0;i<dataLen;i++){
//				_sPacks[j].data[i]
//					= sendX_slave[rNeighbor_slave[kk+i]-cell_slave[0]];
//			}
//			ownNeiSendIdx[_sPacks[j].dst_id]+=dataLen;
//			_sPacks[j].indM = dataLen;
//		}else{
//			for(i=0;i<endIdx-startIdx;i++){
//				_sPacks[j].data[i]
//					= sendX_slave[rNeighbor_slave[kk+i]-cell_slave[0]];
//			}
//			ownNeiSendIdx[_sPacks[j].dst_id]+=endIdx-startIdx;
//			_sPacks[j].indM = endIdx-startIdx;
//		}
//	}
//
//	transform_data();
//	nonZeroNum=0;
//	for(i=0;i<total_send_pcg;i++){
//		if(_rPacks[i].src_id>myId){
//			for(j=0;j<_rPacks[i].indM;j++){
//				recvX_slave[nonZeroNum] = _rPacks[i].data[j];
//				nonZeroNum++;
//			}
//		}
//	}
//
//    //计算稠密块中上三角
//    if(blockLen>0){
//		startIdx = blockStarts_slave[6];
//	    DMA_Get(&sLower_slave[0],&upperReordered[startIdx],
//		  			blockLen*sizeof(swFloat));
//		DMA_Get(&sNeighbor_slave[0],&ownReordered[startIdx],
//		  			blockLen*sizeof(swInt));
//  		for(i=0;i<nonZeroNum;i++){
//  			b_slave[sNeighbor_slave[i]-cell_slave[0]]
//				+= sLower_slave[i]*recvX_slave[i];
//  		}
//	}
//
//	startIdx = blockStarts_slave[2];
//	blockLen = blockStarts_slave[3]-startIdx;
//  	DMA_Get(&diagUpper_slave[0],&upperReordered[startIdx],
//  				blockLen*sizeof(swFloat));
//	DMA_Get(&diagOwner_slave[0],&ownReordered[startIdx],
//			  blockLen*sizeof(swInt));
//  	DMA_Get(&diagNeighbor_slave[0],&neiReordered[startIdx],
//  				blockLen*sizeof(swInt));
//	for(k=0;k<blockStarts_slave[3]-blockStarts_slave[2];k++){
//		b_slave[diagOwner_slave[k]-cell_slave[0]]
//			+= diagUpper_slave[k]
//			*  sendX_slave[diagNeighbor_slave[k]-cell_slave[0]];
//	}
//
//  	DMA_Get(&diagUpper_slave[0],&lowerReordered[startIdx],
//			  blockLen*sizeof(swFloat));
//	for(k=0;k<blockStarts_slave[3]-blockStarts_slave[2];k++){
//		b_slave[diagNeighbor_slave[k]-cell_slave[0]]
//			+= diagUpper_slave[k]
//			*  sendX_slave[diagOwner_slave[k]-cell_slave[0]];
//	}
//
//  	DMA_Get(&recvX_slave[0],&diagReordered[cell_slave[0]],
//			  (cell_slave[1]-cell_slave[0])*sizeof(swFloat));
//	for(k=0;k<cell_slave[1]-cell_slave[0];k++){
//		b_slave[k] += recvX_slave[k]*sendX_slave[k];
//	}
//
//	DMA_Put(&b_MLBSlave[cell_slave[0]],&b_slave[0],
//				(cell_slave[1]-cell_slave[0])*sizeof(swFloat));

}

#include "DS_edge2VertexIter_host.h"
#include <stdlib.h>
#include <stdio.h>

extern "C"
{
#include <athread.h>
void slave_DS_edge2VertexIter_slave(DS_e2VParameters *e2VParas)
}

// Constructor
DS_edge2VertexIter::DS_edge2VertexIter(DSParameters* DSParas)
{
	swInt isXExist,spIndex,dimension,vertexNumber,edgeNumber;
	isXExist     = 1;
	spIndex      = 0;
	dimension    = 1;
	vertexNumber = DSParas->vertexNumber;
	edgeNumber   = DSParas->edgeNumber;
	this->e2VParas = {DSParas,NULL,NULL,0,isXExist,dimension};
	this->segStarts    = DSParas->segStarts;
	this->edgeStarts   = DSParas->edgeStarts;
	this->segNum       = DSParas->segNum;

	initOwnNeiSendList(DSParas);
	// Defined in PROJECT_DIR/RL_MPI/register.C
	initSendList(ownNeiSendList, dataList);
}

// Iteration
void DS_edge2VertexIter::edge2VertexIteration_host(Arrays* upper,
			Arrays* lower, Arrays* diag, Arrays* vertex);
{
	athread_init();
	this->e2VParas->upper  = upper;
	this->e2VParas->lower  = lower;
	this->e2VParas->diag   = diag;
	this->e2VParas->vertex = vertex;
	for(int spIndex=0;spIndex<segNum;spIndex++)
	{
		// Defined in PROJECT_DIR/RL_MPI/register.C
		initTable(spIndex);
		this->e2VParas.spIndex = spIndex;
		__real_athread_spawn((void*)DS_edge2VertexIter_slave,
					&this->e2VParas);
		// Computing in main core
		for(int iseg=spIndex*BLOCKNUM64K;iseg<(spIndex+1)*BLOCKNUM64K;
					iseg++)
		{
			// Get the start neighbor index of sparser(off-diagonal) block
			int minCol = this->segStarts[(spIndex+1)*BLOCKNUM64K];
			swInt dimension = vertex->fArrayDims[0];
			swInt vertexNumber = vertex->fArrayNum;
			swInt edgeNumber   = vertex->fArrayNum;
        	swFloat* x = vertex->floatArrays[1];
        	swFloat* b      = vertex->floatArrays[0];
        	swFloat* lower  = lower->floatArrays[0];
        	swFloat* upper  = upper->floatArrays[1];
        	swInt* owner    = this->e2VParas->DSParas->owner;
        	swInt* neighbor = this->e2VParas->DSParas->neighbor;

			for(int iedge=this->segStarts[iseg];
						this->iedge<segStarts[iseg+1];iedge++)
			{
				if(neighbor[iedge]>=minCol)
				{
					b[owner[iedge]]    += upper[iedge]*x[neighbor[iedge]];
					b[neighbor[iedge]] += lower[iedge]*x[owner[iedge]];
				}
			}
		}
		athread_join();
		// Defined in PROJECT_DIR/RL_MPI/register.C
		destroyTable(spIndex);
	}
	athread_halt();
}

// Init the send/receive topology during the register communication
void DS_edge2VertexIter::initOwnNeiSendList(DSParas)
{
	swInt maxEdges,segNum;
	maxEdges  = DSParas->maxEdges;
	segNum    = DSParas->segNum;

	swInt *owner,*neighbor,*segStarts,*edgeStarts;
	owner      = DSParas->owner;
	neighbor   = DSParas->neighbor;
	segStarts  = DSParas->segStarts;
	edgeStarts = DSParas->edgeStarts;
	
	dataList.resize(segNum*BLOCKNUM64K);
	ownNeiSendList
		= (int*)malloc(segNum*BLOCKNUM64K*BLOCKNUM64K*sizeof(int));
	for(int i=0;i<segNum*BLOCKNUM64K*BLOCKNUM64K;i++){ownNeiSendList[i]=0;}

	// Get the count of edges of every column block.
	for(int iseg=0;iseg<segNum*BLOCKNUM64K;iseg++)
	{
		int blockIdx = iseg/BLOCKNUM64K;
		for(int iedge=edgeStarts[iseg];iedge<edgeStarts[iseg+1];iedge++)
		{
			for(int col=blockIdx*BLOCKNUM64K;
						col<(blockIndex+1)*BLOCKNUM64K;col++)
			{
				if(neighbor[iedge]<segStarts[col+1]
							&& neighbor[iedge>=segStarts[col]])
				{
					int idx = iseg*BLOCKNUM64K+col-blockIdx*BLOCKNUM64K;
					ownNeiSendList[idx]++;
				}
			}
		}
	}

	// Get the count of packages during the register communication.
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
			int packNum = (ownNeiSendList[idx]-1)/6+1;
			int symmIdx = (col+blockIdx*BLOCKNUM64K)*BLOCKNUM64K+irow;
			ownNeiSendList[idx] = packNum;
			ownNeiSendList[symmIdx] = packNum;
		}
	}

	// Get the count of efficient data in every package.
	for(int iseg=0;iseg<segNum*BLOCKNUM64K;iseg++)
	{
		int blockIdx = iseg/BLOCKNUM64K;
		int irow = iseg-blockIdx*BLOCKNUM64K;
		for(int icol=0;icol<BLOCKNUM64K;icol++)
		{
			int idx = iseg*BLOCKNUM64K+icol;
			int packNum = ownNeiSendList[idx];
			while(packNum>6)
			{
				dataList[iseg].push_back(6);
				packNum -=6;
			}
			dataList[iseg].push_back(packNum);
		}
	}
}



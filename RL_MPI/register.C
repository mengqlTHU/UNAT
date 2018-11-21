#include "register.H"
#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>
#include "swMacro.h"
using namespace std;
extern "C"{
#include "athread.h"
#include <algorithm>
void slave_initTable(Schedule *data);
}
void initSendList(int *dataSendList, int mshBlockNum)
{
	printf("start Init SendList...\n");
	schedule_data = new Schedule[mshBlockNum];;
	for(int blockIdx=0;blockIdx<mshBlockNum;blockIdx++){
		RlmpiInitializer reg;
		vector<vector<ThreadID> > sendList(BLOCKNUM64K);
		for(int i=0;i<BLOCKNUM64K;i++){
			int row = blockIdx*BLOCKNUM64K+i;
			for(int j=0;j<BLOCKNUM64K;j++){
				int index = row*BLOCKNUM64K+j;
				for(int k=0;k<dataSendList[index];k++){
					sendList[i].push_back(j);
				}
			}
		}
		reg.init(sendList);
		reg.copyinfo(&schedule_data[blockIdx]);
	}

	printf("finish Init SendList...\n");
}

void initTable(int index){
	schedule_data[index].destroy = 0;
	__real_athread_spawn((void*)slave_initTable,&schedule_data[index]);
	athread_join();
}

void destroyTable(int index){
	schedule_data[index].destroy = 1;
	__real_athread_spawn((void*)slave_initTable,&schedule_data[index]);
	athread_join();
}

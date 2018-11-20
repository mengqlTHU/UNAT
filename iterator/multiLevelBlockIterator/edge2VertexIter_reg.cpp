#include "edge2VertexIter_reg.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>
using namespace std;
extern "C"{
#include "athread.h"
#include <algorithm>
void slave_initTable(Schedule *data);
}
void edge2VertexIteration_reg(int *dataSendList, int mshBlockNum)
{
	printf("start Init SendList...\n");
	schedule_data = new Schedule[mshBlockNum];;
	for(int blockIdx=0;blockIdx<mshBlockNum;blockIdx++){
		RlmpiInitializer reg;
		vector<vector<ThreadID> > sendList(SPNUMs);
		for(int i=0;i<SPNUMs;i++){
			int row = blockIdx*SPNUMs+i;
			for(int j=0;j<SPNUMs;j++){
				int index = row*SPNUMs+j;
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

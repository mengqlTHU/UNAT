#ifndef REGISTER_H
#define REGISTER_H

#include "RlmpiInitializer.hxx"
#ifdef __cplusplus
extern "C"{
#endif
void initSendList(int *dataSendList, int mshBlockNum);
void initTable(int index);
void destroyTable(int index);
Schedule *schedule_data;
#ifdef __cplusplus
}
#endif
#endif

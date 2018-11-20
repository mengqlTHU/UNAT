#ifndef EDGE2VERTEXITER_H
#define EDGE2VERTEXITER_H

#ifdef __cplusplus
extern "C"{
#endif
void edge2VertexIteration_reg(int *dataSendList, int mshBlockNum);
void initTable(int index);
void destroyTable(int index);
Schedule *schedule_data;
#ifdef __cplusplus
}
#endif
#endif

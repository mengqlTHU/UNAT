#include <slave.h>
#include <dma.h>

#ifndef SLAVEUTILS_H
#define SLAVEUTILS_H

#define ALIGNED(addr) ((((unsigned long)((void*)addr-1)>>5)+1)<<5)

#define INIT_LDM_SPACE(ldm_size_) \
; \
int ldm_size = ldm_size_; \
char ldm_space[ldm_size]; \
char* ldm_space_end = ldm_space;

#define LDM_NEW( ptr, type, length ) \
{ \
	ptr = ALIGNED(ldm_space_end); \
	if((char*)ptr - ldm_space \
				+ sizeof(type)*length >= ldm_size ) \
	{ \
		printf("exceed LDM space! \n"); \
		exit(-1); \
	} \
	else \
	{ \
		ldm_space_end = (char*)ptr + sizeof(type)*length; \
	} \
}

#define A_DMA_GET_SET(da,mode,len,re_addr) \
{ \
    dma_set_op(&da,DMA_GET); \
    dma_set_mode(&da,mode); \
    dma_set_size(&da,len); \
    dma_set_reply(&da,re_addr); \
}

#define A_DMA_GET_RUN(da,src,dest) \
{ \
    dma(da,src,dest); \
}

#define A_DMA_PUT_SET(da,mode,len,re_addr) \
{ \
    dma_set_op(&da,DMA_PUT); \
    dma_set_mode(&da,mode); \
    dma_set_size(&da,len); \
    dma_set_reply(&da,re_addr); \
}

#define A_DMA_PUT_RUN(da,src,dest) \
{ \
    dma(da,dest,src); \
}


typedef unsigned long DMA_Status;

void DMA_Get(void* dest, const void* source, const int size);

void DMA_IGet(void* dest, const void* source, const int size, 
        DMA_Status* status);

void DMA_Put(void* dest, const void* source, const int size);

void DMA_IPut(void* dest, const void* source, const int size, 
        DMA_Status* status);

void DMA_Wait( DMA_Status* status, const int num);

void athread_wait(int* reply, int count);

#endif

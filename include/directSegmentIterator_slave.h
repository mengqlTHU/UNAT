#ifndef DIRECTSEGMENTITERATOR_SLAVE_H
#define DIRECTSEGMENTITERATOR_SLAVE_H

#define COPYARRAY(slave, host, length, startIdx) \
{ \
	fArrayNum    = host->fArrayNum; \
	slave->fArrayNum   = fArrayNum; \
	slave->fArraySizes = length; \
	LDM_NEW(slave->fArrayDims, swInt, fArrayNum); \
	LDM_NEW(slave->fArrayInOut, swInt, fArrayNum); \
	DMA_Get(slave->fArrayDims,host->fArrayDims,fArrayNum*sizeof(swInt)); \
	DMA_Get(slave->fArrayInOut,host->fArrayInOut,fArrayNum*sizeof(swInt));\
}
//	slave->floatArrays = (swFloat**)ldm_malloc(fArrayNum*sizeof(swFloat*));\
//	for(i=0;i<fArrayNum;i++) \
//	{ \
//		slave->floatArrays[i] = (swFloat*)ldm_malloc \
//			(length*slave->fArrayDims[i]*sizeof(swFloat)); \
//		DMA_Get(slave->floatArrays[i], \
//			&host->floatArrays[i][slave->fArrayDims[i]*startIdx], \
//			length*slave->fArrayDims[i]*sizeof(swFloat)); \
//	}  \
}

#endif


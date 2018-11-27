#ifndef SWMACRO_H
#define SWMACRO_H

//#include "slaveUtils.h"

#define BLOCKNUM64K 64
#define EPS 1e-6
#define DEBUG

#ifdef DEBUG
#define LOG(format,...) printf("File: "__FILE__",Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)
#else
#define LOG(format,...)
#endif


typedef int swInt;
typedef int swInt32;
typedef long swInt64;
typedef double swFloat;
typedef float swFloat32;
typedef double swFloat64;

#endif

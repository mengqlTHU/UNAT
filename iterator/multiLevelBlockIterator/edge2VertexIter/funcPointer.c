#include "funcPointer.h"
#include <stdlib.h>
#include <stdio.h>
#include <athread.h>


extern SLAVE_FUN(spMV)(MLBFunParameters *MLBFunParas);

PF funcPointer()
{
	return slave_spMV;
}

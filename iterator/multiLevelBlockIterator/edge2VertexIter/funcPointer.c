#include "funcPointer.h"
#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include "userFunc_host.h"


extern SLAVE_FUN(spMV)(MLBFunParameters *MLBFunParas);

PF funcPointer_host()
{
	return spMV;
}

PF funcPointer_slave()
{
	return slave_spMV;
}

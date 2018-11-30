#include "funcPointer.h"
#include <stdlib.h>
#include <stdio.h>
#include <athread.h>
#include "userFunc_host.h"
#include "userFuncUnsymm_host.h"


extern SLAVE_FUN(spMV)(MLBFunParameters *MLBFunParas);
extern SLAVE_FUN(spMVUnsymm)(MLBFunParameters *MLBFunParas);
extern SLAVE_FUN(integrate)(MLBFunParameters *MLBFunParas);
extern SLAVE_FUN(integrateUnsymm)(MLBFunParameters *MLBFunParas);

PF funcPointer_host(swInt flag)
{
	switch (flag)
	{
		case 0:
			return spMVUnsymm;
		case 1:
			return integrateUnsymm;
		default:
			return NULL;
	}
}

PF funcPointer_slave(swInt flag)
{
	switch (flag)
	{
		case 0:
			return slave_spMVUnsymm;
		case 1:
			return slave_integrateUnsymm;
		default:
			return NULL;
	}
}

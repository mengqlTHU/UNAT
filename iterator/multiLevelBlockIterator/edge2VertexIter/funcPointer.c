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
			return spMV;
		case 1:
			return integrate;
		case 2:
			return spMVUnsymm;
		case 3:
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
			return slave_spMV;
		case 1:
			return slave_integrate;
		case 2:
			return slave_spMVUnsymm;
		case 3:
			return slave_integrateUnsymm;
		default:
			return NULL;
	}
}

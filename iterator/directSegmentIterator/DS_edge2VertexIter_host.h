#ifndef DS_EDGE2VERTEXITER_H
#define DS_EDGE2VERTEXITER_H

#include <vertor>
#include "swMacro.h"
#include "iterator_struc.h"
#include "DS_host.h"

// ************************************************************************
//                      class DS_edge2VertexIter
// ************************************************************************

class DS_edge2VertexIter
{
	private:
		// send/receive topology during the register communication
		// dimension(1): segNum*BLOCKNUM64K
		// dimension(2): BLOCKNUM64K
		swInt* ownNeiSendList;
		// the efficient data count of every package
		// dimension(1): segNum*BLOCKNUM64K
		// dimension(2): depends on the dimension(2) of ownNeiSendList
		vector<vector<int> > dataList;
		swInt* segStarts;
		swInt* edgeStarts;
		swInt  segNum;
		DS_e2VParameters e2VParas;

		void initOwnNeiSendList(DSParameters* DSParas);

	public:
		DS_edge2VertexIter(DSParameters* DSParas);
		~DS_edge2VertexIter(){};

		void edge2VertexIteration_host(Arrays* upper,Arrays* lower,
					Arrays* diag, Arrays* vertex);
}


#endif

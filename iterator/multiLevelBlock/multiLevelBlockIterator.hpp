#ifndef UNAT_MULTILEVELBLOCKITERATOR_HPP
#define UNAT_MULTILEVELBLOCKITERATOR_HPP

#include <iostream>
#include "swMacro.h"
#include "iterator.hpp"
#include "iterator.h"
#include "RlmpiInitializer.hxx"
//#include "topology.H"
//using namespace std;
namespace UNAT
{

// ------------------------------------------------------------------------
//                     class MultiLevelBlockIterator
// ------------------------------------------------------------------------

class MultiLevelBlockIterator : public Iterator
{
	private:
		swInt  _cpeBlockNum;
		swInt  _mshBlockNum;
		swInt  _mtxBlockNum;
		swInt* _blockStarts;
		swInt* _blockStartsUnsymm;
		swInt* _vertexStarts;
		swInt  _maxXNum;
		swInt  _maxCells;
		swInt  _maxEdges;
		swInt  _maxEdgesUnsymm;
		swInt* _owner;
		swInt* _neighbor;
		swInt* _postEdgeOrder;
		swInt* _postVertexOrder;
		// 此处保存一份未重排的COO形式拓扑拷贝
		swInt* _firstEdgeVertices;
		swInt* _vertexNeighbours;
		swInt* _recvStarts;

		// Register communication topology
		RlmpiInfo* _schedule_data;

		void reorderVertexArray(swFloat* array);
		void reorderEdgeArrayUnsymm(swFloat* array);
		void MLBReorder(Topology &topo, swInt ref);

		void initOwnNeiSendList();
		void writeTopology();
	public:
		// Constructors
		MultiLevelBlockIterator(Topology &topo, swInt *vertexWeights,
					swInt* edgeWeights);
		// Deconstructors
		~MultiLevelBlockIterator(){};

		virtual void reorderEdges(swInt* startVertices, swInt* endVertices,
					swInt edgeNumber, swInt vertexNumber);
		virtual void reorderNeighbors(swInt* firstEdgeVertices,
					swInt* vertexNeighbours,
					swInt edgeNumber, swInt vertexNumber);
		virtual void edge2VertexIteration(Arrays* backEdgeData,
					Arrays* frontEdgeData, Arrays* selfConnData,
					Arrays* vertexData,
					e2v_hostFunPtr fun_host, e2v_slaveFunPtr fun_slave);
		virtual void vertex2EdgeIteration(Arrays* neighbourData,
					Arrays* vertexData,
					v2e_hostFunPtr, v2e_slaveFunPtr);
		virtual void reorderEdgeData(Arrays* backEdgeData,
					Arrays* frontEdgeData);
		void reorderEdgeDataUnsymm(Arrays* edgeData);
		virtual void reorderVertexData(Arrays* edgeData);

		swInt getCpeBlockNum(){return this->_cpeBlockNum;};
		swInt getMshBlockNum(){return this->_mshBlockNum;};
		swInt getMtxBlockNum(){return this->_mtxBlockNum;};
		swInt getMaxXNum(){return this->_maxXNum;};
		swInt getMaxCells(){return this->_maxCells;};
		swInt getMaxEdges(){return this->_maxEdges;};
		swInt getMaxEdgesUnsymm(){return this->_maxEdgesUnsymm;};
		swInt* getBlockStarts(){return this->_blockStarts;};
		swInt* getBlockStartsUnsymm(){return this->_blockStartsUnsymm;};
		swInt* getVertexStarts(){return this->_vertexStarts;};
};

}

#endif

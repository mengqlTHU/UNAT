#ifndef MULTILEVELBLOCKITERATOR_H
#define MULTILEVELBLOCKITERATOR_H

#include <iostream>
#include "swMacro.h"
#include "iterator.H"
#include "topology.H"

// ------------------------------------------------------------------------
//                     class MultiLevelBlockIterator
// ------------------------------------------------------------------------

class MultiLevelBlockIterator : public Iterator
{
	private:
		swInt    _cpeBlockNum;
		swInt    _mshBlockNum;
		swInt    _mtxBlockNum;
		swInt*   _blockStarts;
		swInt*   _vertexStarts;
		swInt*   _owner;
		swInt*   _neighbor;
		swInt    _maxXNum;
		swInt    _maxCells;
		swInt    _maxEdges;

		void reorderVertexArray(swFloat* array);
	public:
		// Constructors
		MultiLevelBlockIterator(){};
		MultiLevelBlockIterator(Topology &topo) : Iterator(topo) {
			cout<<"MultiLevelBlockIterator"<<endl;
		};
		// Deconstructors
		~MultiLevelBlockIterator(){};

		void reorderEdgesFromEdge(swInt* startVertices, swInt* endVertices,
					swInt edgeNumber, swInt vertexNumber);
		void reorderEdgesFromVertex(swInt* accuVertexEdgeNumbers,
					swInt* vertexNeighbours,
					swInt edgeNumber, swInt vertexNumber);
		void edge2VertexIteration(Arrays* edgeData,
					Arrays* vertexData,
					void (*operatorFunPointer)
					(Arrays*, Arrays*, swInt*, swInt*));
		void vertex2EdgeIteration(Arrays* neighbourData,
					Arrays* vertexData,
					void (*operatorFunPointer)
					(Arrays*, Arrays*, Topology *topo));
		void reorderEdgeData(Arrays* edgeData);
		void reorderVertexData(Arrays* edgeData);

		swInt getCpeBlockNum(){return this->_cpeBlockNum;};
		swInt getMshBlockNum(){return this->_mshBlockNum;};
		swInt getMtxBlockNum(){return this->_mtxBlockNum;};
		swInt getMaxXNum(){return this->_maxXNum;};
		swInt getMaxCells(){return this->_maxCells;};
		swInt getMaxEdges(){return this->_maxEdges;};
		swInt* getBlockStarts(){return this->_blockStarts;};
		swInt* getVertexStarts(){return this->_vertexStarts;};
		swInt* getOwner(){return this->_owner;};
		swInt* getNeighbor(){return this->_neighbor;};
};

#include "multiLevelBlockIterator.C"

#endif

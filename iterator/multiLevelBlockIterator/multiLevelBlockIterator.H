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
		swInt    _maxXNum;
		swInt    _maxCells;
		swInt    _maxEdges;

		void reorderVertexArray(swFloat* array);
	public:
		// Constructors
		MultiLevelBlockIterator(Topology &topo);
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
					(MLBFunParameters *MLBFunParas));
		void vertex2EdgeIteration(Arrays* neighbourData,
					Arrays* vertexData,
					void (*operatorFunPointer)
					(MLBFunParameters *MLBFunParas));
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
};

#include "multiLevelBlockIterator.C"

#endif

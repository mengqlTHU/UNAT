#ifndef ITERATOR_H
#define ITERATOR_H

#include "iterator_struct.h"
#include <stdlib.h>
#include <iostream>
#include <map>
#include "swMacro.h"
#include "topology.H"

// ------------------------------------------------------------------------
//                             class iterator
// ------------------------------------------------------------------------
class Iterator
{
	private:
		Topology* _topo;
		map<swInt, swInt> _edgeMap;
		map<swInt, swInt> _vertexMap;

	public:

		// Constructors
		Iterator(){};
		Iterator(Topology &topo)
		{
			cout<<"Iterator constructor"<<endl;
			this->_topo = topo.clone();
		};

		// Deconstructors
		~Iterator(){};

		void reformInnerTopology()
		{
			reorderEdgesFromEdge(this->_topo->getStartVertices(),
						this->_topo->getEndVertices(),
						this->_topo->getEdgeNumber(),
						this->_topo->getVertexNumber());
			reorderEdgesFromVertex(this->_topo->getAccuVertexEdgeNumbers(),
						this->_topo->getVertexNeighbours(),
						this->_topo->getEdgeNumber(),
						this->_topo->getVertexNumber());
		};

		Topology* getTopology(){return this->_topo;};
		map<swInt, swInt>& getEdgeMap(){return this->_edgeMap;};
		map<swInt, swInt>& getVertexMap(){return this->_vertexMap;};

		virtual void reorderEdgesFromEdge(swInt* startVertices,
					swInt* endVertices,
					swInt edgeNumber, swInt vertexNumber) = 0;
		virtual void reorderEdgesFromVertex(swInt* accuVertexEdgeNumbers,
					swInt* vertexNeighbours,
					swInt edgeNumber, swInt vertexNumber) = 0;
		virtual void reorderEdgeData(Arrays* edgeData) = 0;
		virtual void reorderVertexData(Arrays* edgeData) = 0;

		virtual void edge2VertexIteration(Arrays* edgeData,
					Arrays* vertexData,
					void (*operatorFunPointer)
					(MLBFunParameters *MLBFunParas)) = 0;
		virtual void vertex2EdgeIteration(Arrays* neighbourData,
					Arrays* vertexData,
					void (*operatorFunPointer)
					(MLBFunParameters *MLBFunParas)) = 0;
};

#endif

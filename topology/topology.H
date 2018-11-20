#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <stdlib.h>

using namespace std;
#include "swMacro.h"

// ------------------------------------------------------------------------
//                         class Topology
// ------------------------------------------------------------------------
class Topology{
	private:
		swInt  _vertexNumber;
		swInt  _edgeNumber;
		swInt* _startVertices;
		swInt* _endVertices;
		swInt* _startVertexNumbers;
		swInt* _accuStartVertexNumbers;
		swInt* _firstEdgeVertices;
		swInt* _vertexNeighbours;
		swInt* _vertexEdgeNumbers;
		swInt* _accuVertexEdgeNumbers;
		void EdgeBasedInit();
		void VertexBasedInit();
		void copy(const Topology &topo);

	public:
		// Constructors
		Topology();
		// Copy constructors
		Topology(const Topology &topo);
		static Topology constructFromEdge(swInt *startVertices,
					swInt *endVertices, swInt edgeNumber, bool copy=false);
		static Topology constructFromVertex(swInt *accuVertexEdgeNumbers,
					swInt *vertexNeighbours,swInt vertexNumber,
					bool copy=false);
		
		// Deconstructors
		~Topology();

		swInt  getVertexNumber();
		swInt  getEdgeNumber();
		swInt* getStartVertices();
		swInt* getEndVertices();
		swInt* getStartVertexNumbers();
		swInt* getAccuStartVertexNumbers();
		swInt* getFirstEdgeVertices();
		swInt* getVertexNeighbours();
		swInt* getVertexEdgeNumbers();
		swInt* getAccuVertexEdgeNumbers();

		Topology& operator=(const Topology& topo);
		Topology* clone() const {return new Topology(*this);};

		// TODO: 对于未按顺序排列的LDU方式存储拓扑，提供排序接口，
		// 返回类型为LDUMatrix，CSR同理
		// 提供LDU拓扑与CSR拓扑之间的相互转换接口
		void addEdge();
		void addVertex();
		void removeEdge();
		void removeVertex();
		void transpose();

		void sortAndCompress();
		void edgeBasedToVertexBased();
		void vertexBasedToEdgeBased();
};

#include "topology.C"

#endif

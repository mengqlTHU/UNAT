/*#########################################################
################## HERE IS THE HEADER ######################
Direct segment describes a decomposing method that decompose
the unstructure data set directly with the original
sequence the data in. 
Two level irregular tile was used implicitly to meet load
balance and efficient cache usage.
Data swapping between parallel threads was used to simulate
shared cache.
Minimize the topology matrix bandwidth or block the segments 
with high cohesion to ensure high efficiency of the Iterator 
##########################################################*/
#ifndef UNAT_DIRECTSEGMENTITERATOR_HPP
#define UNAT_DIRECTSEGMENTITERATOR_HPP

#include "iterator.hpp"
#include "RlmpiInitializer.hxx"

namespace UNAT
{

class DirectSegmentIterator
	:
	public Iterator
{
	// Data to describe segmentation
	swInt segNum_;
	swInt subSegNum_;
	swInt maxRowEdges_;
	swInt *segStarts_;
	swInt *edgeStarts_;
	swInt *accuColNum_;
	swInt *colSegNum_;
	swInt *recvStarts_;
	swInt *wrtStarts_;

	// Data for computing in master core
	swInt  *rOwner_;
	swInt  *rNeighbor_;
	Arrays *rBackEdgeData_;
	Arrays *rFrontEdgeData_;
	Arrays *rVertexData_;

	// topology of segment connectivity
	Topology* segConnetion_;
	// egde Number between segment
	swInt* segEdgeNum_;
	// edge neighbor segment index
	swInt* edgeNeiSeg_;

	// Data for storing the topology of slave cores in register comm.
	RlmpiInfo* schedule_data_;

    // send/receive topology during the register communication
	// dimension(1): segNum*BLOCKNUM64K
	// dimension(2): BLOCKNUM64K
	swInt* ownNeiSendList;
	// the efficient data count of every package
	// dimension(1): segNum*BLOCKNUM64K
	// dimension(2): depends on the dimension(2) of ownNeiSendList
    vector<vector<int> > dataList;

    // Array two level segmentation
	void decomposeArray(swInt subSegNum, swInt subSegDataLimit,
				swInt* dataWeights, swInt size);	
	void initOwnNeiSendList();
	void initSendList(int *dataSendList, vector<vector<int> > dataList,
				int mshBlockNum);
public:
	// Dummy Constructor
	//DirectSegmentIterator();
	// Constructor and assignment
	DirectSegmentIterator(Topology &topo, swInt* vertexWeights, swInt* edgeWeights);
	// Copy constructors
	DirectSegmentIterator(const DirectSegmentIterator &);
	DirectSegmentIterator& operator=(const DirectSegmentIterator& );
	DirectSegmentIterator* clone() const {return new DirectSegmentIterator(*this);};

	// Deconstructor
	virtual ~DirectSegmentIterator();

	// Iteration processor
	virtual void edge2VertexIteration(Arrays* backEdgeData,
				Arrays* frontEdgeData, Arrays* selfConnData, Arrays* vertexData, 
				e2v_hostFunPtr fun_host, e2v_slaveFunPtr fun_slave);

	virtual void vertex2EdgeIteration(Arrays* neighbourData,
				Arrays* vertexData, v2e_hostFunPtr fun_host, v2e_slaveFunPtr fun_slave);


	// Data to describe segmentation
	swInt getSegNum()
	{ return segNum_; }
	swInt getSubSegNum()
	{ return subSegNum_; }
	swInt* getSegStarts()
	{ return segStarts_; }
	swInt* getEdgeStarts()
	{ return edgeStarts_; }

	// topology of segment connectivity
	Topology* getSegConnection()
	{ return segConnetion_; }
	// egde Number between segment
	swInt* getSegEdgnNum()
	{ return segEdgeNum_; }
	// edge neighbor segment index
	swInt* getEdgeNeiSeg()
	{ return edgeNeiSeg_; }

}; // class DirectSegmentIterator

} // namespace UNAT

#endif // UNAT_DIRECTSEGMENT_HPP

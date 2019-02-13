#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <vector>
#include <cmath>
#include <assert.h>
#include "swMacro.h"
#include "topology.H"
#include "iterator.hpp"
#include "iterator.h"
#include "directSegmentIterator.hpp"
#include "wrappedInterface.h"
#include "athread_switch.h"

using namespace UNAT;

#define NONZERONUM 616320
#define DIMENSION 1

int* readFile(char* name);
void debug(Topology topo);
//void checkResult(swFloat* array1, swFloat* array2, swInt count);
void checkDSI(DirectSegmentIterator& iterator);

double time1, time2, time3, time4;

//extern "C"
//{
//	//define 2 function pointers
//	define_e2v_hostFunPtr(spMV_host)
//	{
//		//selfConn computation
//		swFloat* diag	= accessArray(selfConnData, 0);
//		swFloat* x		= accessArray(vertexData, 0);
//		swFloat* b		= accessArray(vertexData, 1);
//		
//		swInt vertexNum = getArraySize(selfConnData);
//		swInt ivertex;
//		for( ivertex = 0; ivertex < vertexNum; ivertex++)
//		{
//			b[ivertex] = diag[ivertex]*x[ivertex];
//		}
//		
//		//frontEdge computation
//		swFloat* upper	= accessArray(frontEdgeData, 0);
//		swInt edgeNumber = getArraySize( frontEdgeData );
//		swInt iedge;
//		for( iedge = 0; iedge < edgeNumber; iedge++)
//			b[startVertices[iedge]] += upper[iedge]*x[endVertices[iedge]];
//
//		//backEdge computation
//		swFloat* lower	= accessArray(backEdgeData, 0);
//		edgeNumber = getArraySize( backEdgeData );
//		for( iedge = 0; iedge < edgeNumber; iedge++)
//			b[endVertices[iedge]] += lower[iedge]*x[startVertices[iedge]];
//	}
//
//	define_e2v_slaveFunPtr(spMV_slave)
//	{
//		//selfConn computation
//		swFloat* diag	= accessArray(selfConnData, 0);
//		swFloat* x		= accessArray(vertexData, 0);
//		swFloat* b		= accessArray(vertexData, 1);
//		
//		swInt vertexNum = getArraySize(selfConnData);
//		swInt ivertex;
//		for( ivertex = 0; ivertex < vertexNum; ivertex++)
//		{
//			b[ivertex] = diag[ivertex]*x[ivertex];
//		}
//		
//		//frontEdge computation
//		swFloat* upper	= accessArray(frontEdgeData, 0);
//		swInt edgeNumber = getArraySize(frontEdgeData);
//		swInt iedge;
//		for( iedge = 0; iedge < edgeNumber; iedge++)
//			b[startVertices[iedge]] += upper[iedge]*x[endVertices[iedge]];
//
//		//backEdge computation
//		swFloat* lower	= accessArray(backEdgeData, 0);
//		edgeNumber = getArraySize( backEdgeData );
//		for( iedge = 0; iedge < edgeNumber; iedge++)
//			b[endVertices[iedge]] += lower[iedge]*x[startVertices[iedge]];
//	}
//}

int main()
{
	char owner[] = "owner_of";
	char neighbor[] = "neighbour_of";
	swInt *rowAddr = readFile(owner);
	swInt *colAddr = readFile(neighbor);
	Topology* topo = Topology::constructFromEdge(rowAddr,colAddr,NONZERONUM);
	int dims = DIMENSION;

	swFloat* lower
		= (swFloat*)malloc(sizeof(swFloat)*topo->getEdgeNumber()*dims);
	swFloat* lower_DSI
		= (swFloat*)malloc(sizeof(swFloat)*topo->getEdgeNumber()*dims);
	swFloat* upper
		= (swFloat*)malloc(sizeof(swFloat)*topo->getEdgeNumber()*dims);
	swFloat* upper_DSI
		= (swFloat*)malloc(sizeof(swFloat)*topo->getEdgeNumber()*dims);
	swFloat* diag
		= (swFloat*)malloc(sizeof(swFloat)*topo->getVertexNumber()*dims);
	swFloat* diag_DSI
		= (swFloat*)malloc(sizeof(swFloat)*topo->getVertexNumber()*dims);
	swFloat* x
		= (swFloat*)malloc(sizeof(swFloat)*topo->getVertexNumber()*dims);
	swFloat* x_DSI
		= (swFloat*)malloc(sizeof(swFloat)*topo->getVertexNumber()*dims);
	swFloat* b
		= (swFloat*)malloc(sizeof(swFloat)*topo->getVertexNumber()*dims);
	swFloat* b_DSI
		= (swFloat*)malloc(sizeof(swFloat)*topo->getVertexNumber()*dims);
	//swFloat* data  = (swFloat*)malloc( sizeof(swFloat)
	//								   *( topo->getEdgeNumber()*2 
	//									  + topo->getVertexNumber() )
	//								 );

	// weights
	std::vector<swInt> cellWeights(topo->getVertexNumber(), 2*dims);
	std::vector<swInt> edgeWeights(topo->getEdgeNumber(), dims+1);


	DirectSegmentIterator iterator(*topo, &cellWeights[0], &edgeWeights[0]);


//	checkDSI(iterator);

	// assign value
	printf("edgeNum: %d\n",topo->getEdgeNumber());
	for(int i=0;i<topo->getEdgeNumber();i++)
	{
		for(int iDim=0;iDim<dims;iDim++)
		{
			lower[i*dims+iDim]
				= ((swFloat)rowAddr[i]+1)/((swFloat)colAddr[i]+1);
			lower_DSI[i*dims+iDim] = lower[i*dims+iDim];
			upper[i*dims+iDim]
				= ((swFloat)colAddr[i]+1)/((swFloat)rowAddr[i]+1);
			upper_DSI[i*dims+iDim] = upper[i*dims+iDim];
		}
	}
	for(int i=0;i<topo->getVertexNumber();i++)
	{
		for(int iDim=0;iDim<dims;iDim++)
		{
			diag[i*dims+iDim]      = i;
			diag_DSI[i*dims+iDim]  = i;
			x[i*dims+iDim]     = (double)(i+1)/(i+2);
			x_DSI[i*dims+iDim] = (double)(i+1)/(i+2);
			b[i*dims+iDim]     = (double)(i+1)/(i+2);
			b_DSI[i*dims+iDim] = (double)(i+1)/(i+2);
		}
	}


	printf("Construct array...\n");
	// calculate with iterator and function pointers
	Arrays backEdgeData, frontEdgeData, 
		   selfConnData, vertexData;
	Arrays refBackEdgeData, refFrontEdgeData, 
		   refSelfConnData, refVertexData;
	// get vertex number
	swInt vertexNum = topo->getVertexNumber();
	swInt edgeNum   = topo->getEdgeNumber();
    constructSingleArray(refBackEdgeData, dims, NONZERONUM, COPYIN,
				lower_DSI);
//	constructEmptyArray(refBackEdgeData);
    constructSingleArray(refFrontEdgeData, dims, NONZERONUM, COPYIN,
				upper_DSI);
//	addSingleArray(refFrontEdgeData, 1, NONZERONUM, COPYIN, lower_DSI);
//    constructSingleArray(refSelfConnData, dims, vertexNum, COPYOUT,
//				diag_DSI);
	constructEmptyArray(refSelfConnData);
//	constructEmptyArray(refVertexData);
    constructSingleArray(refVertexData, dims, vertexNum, COPYOUT, b_DSI);
//    addSingleArray(refVertexData, dims, vertexNum, COPYOUT, b_DSI);

	// single constructor for arrays
//	constructEmptyArray(backEdgeData);
    constructSingleArray(backEdgeData, dims, NONZERONUM, COPYIN, lower);
    constructSingleArray(frontEdgeData, dims, NONZERONUM, COPYIN, upper);
//	addSingleArray(frontEdgeData, 1, NONZERONUM, COPYIN, lower);
//    constructSingleArray(selfConnData, dims, vertexNum, COPYOUT, diag);
//	constructEmptyArray(frontEdgeData);
	constructEmptyArray(selfConnData);
//	constructEmptyArray(vertexData);
    constructSingleArray(vertexData, dims, vertexNum, COPYOUT, b);
//    constructSingleArray(vertexData, dims, vertexNum, COPYIN, x);
//    addSingleArray(vertexData, dims, vertexNum, COPYOUT, b);
//    constructSingleArray(refVertexData, dims, vertexNum, COPYOUT, b);
//    addSingleArray(refVertexData, dims, vertexNum, COPYOUT, b);

	// reference calculation
	// Timer
	getTime(time1);

	swNegSumDiag(&refBackEdgeData,&refFrontEdgeData,&refSelfConnData,
				&refVertexData,rowAddr,colAddr);
//	for(int i=0;i<topo->getVertexNumber();i++)
//	{
//		for(int iDim=0;iDim<dims;iDim++)
//		{
////if(i*dims+iDim==26604) printf("diag:%f,%f,%f\n",b[i*dims+iDim],diag[i*dims+iDim],x[i*dims+iDim]);
//			b[i*dims+iDim]
//				= diag[i*dims+iDim]*x[i*dims+iDim];
//		}
//	}
//	for(int i=0;i<topo->getEdgeNumber();i++)
//	{
//		for(int iDim=0;iDim<dims;iDim++)
//		{
////if(rowAddr[i]*dims+iDim==26604) printf("upper:%d,%f,%f,%f\n",colAddr[i],b[rowAddr[i]*dims+iDim],upper[i*dims+iDim],x[colAddr[i]*dims+iDim]);
//			b[rowAddr[i]*dims+iDim]
//				+= upper[i*dims+iDim]*x[colAddr[i]*dims+iDim];
////if(colAddr[i]*dims+iDim==26604) printf("lower:%d,%f,%f,%f\n",rowAddr[i],b[colAddr[i]*dims+iDim],lower[i*dims+iDim],x[rowAddr[i]*dims+iDim]);
//			b[colAddr[i]*dims+iDim]
//				+= lower[i*dims+iDim]*x[rowAddr[i]*dims+iDim];
//		}
//	}
	// Timer
	getTime(time2);
	printf("CPU Processor Time: %f us\n", (time2-time1)*1000000); 


	CG_init();
	// compute results with iterator
	getTime(time3);
	iterator.edge2VertexIteration( &backEdgeData, &frontEdgeData, 
				&selfConnData, &vertexData,
				swNegSumDiag, slave_swNegSumDiag);
	getTime(time4);
	printf("Slave Core Time: %f us\n", (time4-time3)*1000000); 
	printf("Speed-up: %f\n", (time2-time1)/(time4-time3)); 

	CG_halt();
	checkResult(b_DSI, b, vertexNum*dims);
//	checkResult(diag_DSI, diag, vertexNum*dims);

	free(lower);
	free(upper);
	free(diag);
	free(x);
	free(b);
	free(b_DSI);
	//free(data);
	
	return 0;
}

//void checkResult(swFloat* array1, swFloat* array2, swInt count)
//{
//	cout<<"check result..."<<endl;
//	for(int i=0;i<count;i++)
//	{
//		if(fabs(array1[i]-array2[i])>EPS)
//		{
//			if(fabs(array1[i]-array2[i])/array1[i]>EPS)
//			{
//				cout<<"Error on index["<<i<<"], "
//					<<array1[i]<<", "<<array2[i]<<endl;
//				exit(1);
//			}
//		}
//	}
//	cout<<"The result is correct!"<<endl;
//}

int* readFile(char* name){
	FILE *fp = fopen(name,"r");
	if(fp==NULL){
		return 0;
	}
	char c1[100],c2[100];
	int cellNum;
	while(1)
	{
		fgets(c1,100,fp);
		if(c1[0]=='('){
			cellNum = atoi(c2);
			break;
		}
		fgets(c2,100,fp);
		if(c2[0]=='('){
			cellNum = atoi(c1);
			break;
		}else{
		}
	}
	int *buff = (int*)malloc(sizeof(int)*NONZERONUM);
	int i=0;
	while(i<NONZERONUM)
	{
		fgets(c1,100,fp);
		buff[i] = atoi(c1);
		i++;
	}
	printf("InternalFaceNum = %d\n",i);
	return buff;
}

void debug(Topology topo)
{
	cout<<"edgeNumber: "<<topo.getEdgeNumber()<<endl;
	cout<<"vertexNumber: "<<topo.getVertexNumber()<<endl;

	cout<<"startVertices:"<<endl;
	swInt *startVertices=topo.getStartVertices();
	for(int i=0;i<topo.getEdgeNumber();i++)
	{
		cout<<" "<<startVertices[i];
	}
	cout<<endl;

	cout<<"endVertices:"<<endl;
	swInt *endVertices=topo.getEndVertices();
	for(int i=0;i<topo.getEdgeNumber();i++)
	{
		cout<<" "<<endVertices[i];
	}
	cout<<endl;

	cout<<"startVertexNumbers:"<<endl;
	swInt *startVertexNumbers=topo.getStartVertexNumbers();
	for(int i=0;i<topo.getVertexNumber();i++)
	{
		cout<<" "<<startVertexNumbers[i];
	}
	cout<<endl;

	cout<<"accuStartVertexNumbers:"<<endl;
	swInt *accuStartVertexNumbers=topo.getAccuStartVertexNumbers();
	for(int i=0;i<topo.getVertexNumber()+1;i++)
	{
		cout<<" "<<accuStartVertexNumbers[i];
	}
	cout<<endl;

	cout<<"firstEdgeVertices:"<<endl;
	swInt *firstEdgeVertices=topo.getFirstEdgeVertices();
	for(int i=0;i<topo.getEdgeNumber()*2;i++)
	{
		cout<<" "<<firstEdgeVertices[i];
	}
	cout<<endl;

	cout<<"vertexNeighbours:"<<endl;
	swInt *vertexNeighbours=topo.getVertexNeighbours();
	for(int i=0;i<topo.getEdgeNumber()*2;i++)
	{
		cout<<" "<<vertexNeighbours[i];
	}
	cout<<endl;

	cout<<"vertexEdgeNumbers:"<<endl;
	swInt *vertexEdgeNumbers=topo.getVertexEdgeNumbers();
	for(int i=0;i<topo.getVertexNumber();i++)
	{
		cout<<" "<<vertexEdgeNumbers[i];
	}
	cout<<endl;

	cout<<"accuVertexEdgeNumbers:"<<endl;
	swInt *accuVertexEdgeNumbers=topo.getAccuVertexEdgeNumbers();
	for(int i=0;i<topo.getVertexNumber()+1;i++)
	{
		cout<<" "<<accuVertexEdgeNumbers[i];
	}
	cout<<endl;
}

void checkDSI(DirectSegmentIterator& iterator)
{
	swInt segNum = iterator.getSegNum();
	swInt subSegNum = iterator.getSubSegNum();
	swInt* segStarts = iterator.getSegStarts();
	swInt* edgeStarts = iterator.getEdgeStarts();

	Topology* segConnection = iterator.getSegConnection();
	swInt* segEdgeNum = iterator.getSegEdgnNum();
	swInt* edgeNeiSeg = iterator.getEdgeNeiSeg();

	swInt totalSeg = segNum*subSegNum;

	// print segment info
	printf("##########################################################\n");
	printf("######## Total segment number is %d\n", totalSeg);
	printf("######## Subsegment number is %d\n", subSegNum);
	for(swInt iseg = 0; iseg< segNum; iseg++)
	{
		printf("#### %dth segmet is:\n", iseg);
		printf("vertex starts: ");
		printArray("%d", &segStarts[iseg*subSegNum], subSegNum+1);
		for(swInt isub=0; isub<subSegNum; isub++)
		{
			swInt pSeg = iseg*subSegNum+isub;
			if(segStarts[pSeg]>segStarts[pSeg+1])
			  printf("\n****Error: start %d < end %d\n", segStarts[pSeg], segStarts[pSeg+1]);
		}		
		printf("edge starts: ");
		printArray("%d", &edgeStarts[iseg*subSegNum], subSegNum+1);
		for(swInt isub=0; isub<subSegNum; isub++)
		{
			swInt pSeg = iseg*subSegNum+isub;
			if(edgeStarts[pSeg]>edgeStarts[pSeg+1])
			  printf("\n****Error: start %d < end %d\n", edgeStarts[pSeg], edgeStarts[pSeg+1]);
		}		
	}

	printf("##########################################################\n");

	swInt connetNum = segConnection->getEdgeNumber();
	printf("######## Total segment connetion is %d\n", connetNum);
	swInt* startSeg = segConnection->getStartVertices();
	swInt* endSeg = segConnection->getEndVertices();
	printArray("%d", startSeg, segConnection->getEdgeNumber());
	printArray("%d", endSeg, segConnection->getEdgeNumber());
	printArray("%d", segEdgeNum, segConnection->getEdgeNumber());

	Topology* basicTopo = iterator.getTopology();
	swInt edgeNum = basicTopo->getEdgeNumber();
	swInt* startVertex = basicTopo->getStartVertices();
	swInt* endVertex = basicTopo->getEndVertices();
	

	// calculate segEdgeNum from edgeNeiSeg
	std::vector<swInt> segEdgeNumTmp(totalSeg, 0);
	for(swInt iedge = 0; iedge< edgeNum; iedge++)
	{
		if(edgeNeiSeg[iedge]>=totalSeg || edgeNeiSeg[iedge] < 0)
			printf("****Error: wrong edgeNeiSeg at seg %d\n", edgeNeiSeg[iedge]);
		swInt segRow, segCol;
		for(swInt iseg = 0; iseg<totalSeg; iseg++)
		{
			if(startVertex[iedge] >= segStarts[iseg])
				segRow = iseg;
			if(endVertex[iedge] >= segStarts[iseg])
				segCol = iseg;
		}
		if(segRow != segCol) 
			segEdgeNumTmp[edgeNeiSeg[iedge]]++;
	}
	// check segEdgeNum with segEdgeNumTmp
	std::vector<swInt> segEdgeNumTmp2(totalSeg, 0);
	for(swInt iconn = 0; iconn<connetNum; iconn++)
		segEdgeNumTmp2[endSeg[iconn]] += segEdgeNum[iconn];
	for(swInt iseg = 0; iseg<totalSeg; iseg++)
	{
		if(segEdgeNumTmp2[iseg] != segEdgeNumTmp[iseg])
			printf("****Error: different segEdgeNum for segCol %d: %d to %d\n!", 
						iseg, segEdgeNumTmp2[iseg], segEdgeNumTmp[iseg]);
	}
	// print segment connectivity
	swInt segCount = 0;
	printf("##########################################################\n");
	printf("######## Total segment number is %d\n", segConnection->getVertexNumber());
	if(segConnection->getVertexNumber() != totalSeg )
	{
		printf("****Warning: different segment number!\n");
	}
	for(swInt iseg = 0; iseg< segNum; iseg++)
	{
		printf("-----------------------------\n");
		printf("#### the %dth segment:\n", iseg);
		for(swInt isub=0; isub<subSegNum; isub++)
		{
			swInt pSeg = iseg*subSegNum+isub;
			printf("%dth: ", pSeg);
			while(segCount < segConnection->getEdgeNumber()
						&& startSeg[segCount] == startSeg[segCount+1])
			{
				printf("(%d, %d){%d In%d}, ", 
							startSeg[segCount], 
							endSeg[segCount], 
							segEdgeNum[segCount], 
							segEdgeNumTmp[endSeg[segCount]] );
				segCount++;
			}
			if(segCount < segConnection->getEdgeNumber())
			{
				printf("(%d, %d){%d In%d}, ", 
							startSeg[segCount], 
							endSeg[segCount], 
							segEdgeNum[segCount], 
							segEdgeNumTmp[endSeg[segCount]] );
				segCount++;
			}
			printf("\n");
		}		
	}
	printf("##########################################################\n");
}

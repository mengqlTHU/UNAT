#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include "swMacro.h"
#include "topology.H"
#include "iterator.H"
#include "multiLevelBlockIterator.H"
#include "iterator_struct.h"
#include "funcPointer.h"

#define NONZERONUM 1516800

int* readFile(char* name);
void debug(Topology t);
void (*operatorFunPointer_host)(MLBFunParameters *MLBFunParas);
void (*operatorFunPointer_slave)(MLBFunParameters *MLBFunParas);
//void spMV(MLBFunParameters *MLBFunParas);
void checkResult(swFloat* array1, swFloat* array2, swInt count);
//extern void SLAVE_FUNC(spMV)(MLBFunParameters *MLBFunParas);

int main()
{
	char owner[] = "owner_512000";
	char neighbor[] = "neighbour_512000";
	swInt *rowAddr = readFile(owner);
	swInt *colAddr = readFile(neighbor);
	Topology t = Topology::constructFromEdge(rowAddr,colAddr,NONZERONUM);
//	Topology::constructFromEdge(t,rowAddr,colAddr,NONZERONUM);
//	debug(t);


	swInt dimension = 2;
	swFloat* lower = (swFloat*)malloc
		(sizeof(swFloat)*t.getEdgeNumber()*dimension);
	swFloat* upper = (swFloat*)malloc
		(sizeof(swFloat)*t.getEdgeNumber()*dimension);
	swFloat* diag  = (swFloat*)malloc
		(sizeof(swFloat)*t.getVertexNumber()*dimension);
	swFloat* x     = (swFloat*)malloc
		(sizeof(swFloat)*t.getVertexNumber()*dimension);
	swFloat* b     = (swFloat*)malloc
		(sizeof(swFloat)*t.getVertexNumber()*dimension);
	swFloat* b_mlb = (swFloat*)malloc
		(sizeof(swFloat)*t.getVertexNumber()*dimension);
	swFloat* b_csr = (swFloat*)malloc
		(sizeof(swFloat)*t.getVertexNumber()*dimension);
	swFloat* data  = (swFloat*)malloc
		(sizeof(swFloat)*t.getEdgeNumber()*2*dimension);

	for(int j=0;j<dimension;j++)
	{
		for(int i=0;i<t.getEdgeNumber();i++)
		{
			lower[i+j*t.getEdgeNumber()]
				= (swFloat)(rowAddr[i]+1)/(colAddr[i]+1);
			upper[i+j*t.getEdgeNumber()]
				= (swFloat)(colAddr[i]+1)/(rowAddr[i]+1);
		}
	}
	for(int j=0;j<dimension;j++)
	{
		for(int i=0;i<t.getVertexNumber();i++)
		{
			diag[i+j*t.getVertexNumber()]  = i+1;
			x[i+j*t.getVertexNumber()]     = (double)(i+1)/(i+2);
			b[i+j*t.getVertexNumber()]     = 0;
			b_mlb[i+j*t.getVertexNumber()] = 0;
		}
	}
	for(int j=0;j<dimension;j++)
	{
		for(int i=0;i<t.getEdgeNumber()*2;i++)
		{
			data[i+j*t.getEdgeNumber()*2] = (double)(i+1)/(i+2);
		}
	}

	swInt vertexNumber = t.getVertexNumber();
	swInt edgeNumber   = t.getEdgeNumber();
	struct timeval start, end;
	gettimeofday(&start,NULL);
	for(int j=0;j<dimension;j++)
	{
		for(int i=0;i<vertexNumber;i++)
		{
//		if(i==1)
//		{
//			printf("%d,%f,%f\n",i,diag[i],x[i]);
//		}
			b[i+j*vertexNumber]
				= diag[i+j*vertexNumber]*x[i+j*vertexNumber];
		}
		for(int i=0;i<edgeNumber;i++)
		{
//		if(rowAddr[i]==1) printf("%d,%f,%f\n",i,upper[i],x[colAddr[i]]);
			b[rowAddr[i]+j*vertexNumber]
				+= upper[i+j*edgeNumber] * x[colAddr[i]+j*vertexNumber];
//		if(colAddr[i]==1) printf("%d,%f,%f\n",i,lower[i],x[rowAddr[i]]);
			b[colAddr[i]+j*vertexNumber]
				+= lower[i+j*edgeNumber] * x[rowAddr[i]+j*vertexNumber];
//     	b[rowAddr[i]] += upper[i]*x[colAddr[i]];
//		b[colAddr[i]] += lower[i]*x[rowAddr[i]];
		}
	}
//	gettimeofday(&end,NULL);
//	int timeuse = 1000000*(end.tv_sec-start.tv_sec)
//		+ end.tv_usec-start.tv_usec;
//	printf("CPU Processor Time: %f us\n",(double)timeuse); 

//	swInt* firstEdgeVertices = t.getFirstEdgeVertices();
//	swInt* vertexNeighbor    = t.getVertexNeighbours();
//	for(int i=0;i<t.getEdgeNumber()*2;i++)
//	{
////		if(firstEdgeVertices[i]==6)
////		{
////			printf("%d,%f,%f\n",i,data[i],x[vertexNeighbor[i]]);
////		}
////		b[firstEdgeVertices[i]] += data[i];
//		b[firstEdgeVertices[i]] += data[i]*x[vertexNeighbor[i]];
//	}
	gettimeofday(&end,NULL);
	int timeuse = 1000000*(end.tv_sec-start.tv_sec)
		+ end.tv_usec-start.tv_usec;
	printf("CPU Processor Time: %f us\n",(double)timeuse); 

	operatorFunPointer_host  = funcPointer_host(0);
	operatorFunPointer_slave = funcPointer_slave(0);
	Arrays edgeData
		= {lower, upper, NULL, NULL, t.getEdgeNumber(), dimension};
	Arrays neighbourData
		= { data, NULL, NULL, NULL, t.getEdgeNumber(), dimension};
	Arrays vertexData
		= {b_mlb,    x, diag, NULL, t.getVertexNumber(), dimension};
//	Arrays vertexData    = {b_mlb, NULL, NULL, NULL, t.getVertexNumber()};

	MultiLevelBlockIterator mlbIter(t);
	mlbIter.reformInnerTopology();
	mlbIter.reorderEdgeDataUnsymm(&neighbourData);
	mlbIter.reorderEdgeData(&edgeData);
	mlbIter.reorderVertexData(&vertexData);

	mlbIter.edge2VertexIteration(&edgeData,&vertexData,
				operatorFunPointer_host, operatorFunPointer_slave);
//	mlbIter.vertex2EdgeIteration(&neighbourData,&vertexData,
//				operatorFunPointer_host, operatorFunPointer_slave);

	checkResult(b, vertexData.A4Ptr, t.getVertexNumber()*dimension);

//	swInt accuVertexEdgeNumbers[] = {0,3,6,8,9,12};
//	swInt vertexNeighbours[] = {1,2,4,0,2,4,0,1,4,0,1,3};
//	Topology t2 = Topology::constructFromVertex(accuVertexEdgeNumbers,
//				vertexNeighbours,5);
//	debug(t2);
//
//	Topology t3 = t;
//	t3=t2;
//	debug(t3);

//	Topology* t4=t.clone();
//	debug(*t4);
	return 0;
}

//void spMV(MLBFunParameters *MLBFunParas)
//{
//	cout<<"spMV"<<endl;
////	for(int i=0;i<edgeData->num;i++)
////	{
////		vertexData->A1Ptr[neighbor[i]]
////			+= edgeData->A1Ptr[i]*vertexData->A2Ptr[i];
////	}
//}

void checkResult(swFloat* array1, swFloat* array2, swInt count)
{
	cout<<"check result..."<<endl;
	for(int i=0;i<count;i++)
	{
//		printf("%d,%f,%f\n",i,array1[i],array2[i]);
		if(fabs(array1[i]-array2[i])>EPS)
		{
			if(array2[i]==0)
			{
				if(fabs(array1[i])>1e-14)
				{
					cout<<"Error on index["<<i<<"], "
						<<array1[i]<<", "<<array2[i]<<endl;
					exit(1);
				}
			}else if(fabs(array1[i]-array2[i])/array1[i]>EPS)
			{
				cout<<"Error on index["<<i<<"], "
					<<array1[i]<<", "<<array2[i]<<endl;
				exit(1);
			}
		}
	}
	cout<<"The result is correct!"<<endl;
}

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

void debug(Topology t)
{
	cout<<"edgeNumber: "<<t.getEdgeNumber()<<endl;
	cout<<"vertexNumber: "<<t.getVertexNumber()<<endl;

	cout<<"startVertices:"<<endl;
	swInt *startVertices=t.getStartVertices();
	for(int i=0;i<t.getEdgeNumber();i++)
	{
		cout<<" "<<startVertices[i];
	}
	cout<<endl;

	cout<<"endVertices:"<<endl;
	swInt *endVertices=t.getEndVertices();
	for(int i=0;i<t.getEdgeNumber();i++)
	{
		cout<<" "<<endVertices[i];
	}
	cout<<endl;

	cout<<"startVertexNumbers:"<<endl;
	swInt *startVertexNumbers=t.getStartVertexNumbers();
	for(int i=0;i<t.getVertexNumber();i++)
	{
		cout<<" "<<startVertexNumbers[i];
	}
	cout<<endl;

	cout<<"accuStartVertexNumbers:"<<endl;
	swInt *accuStartVertexNumbers=t.getAccuStartVertexNumbers();
	for(int i=0;i<t.getVertexNumber()+1;i++)
	{
		cout<<" "<<accuStartVertexNumbers[i];
	}
	cout<<endl;

	cout<<"firstEdgeVertices:"<<endl;
	swInt *firstEdgeVertices=t.getFirstEdgeVertices();
	for(int i=0;i<t.getEdgeNumber()*2;i++)
	{
		cout<<" "<<firstEdgeVertices[i];
	}
	cout<<endl;

	cout<<"vertexNeighbours:"<<endl;
	swInt *vertexNeighbours=t.getVertexNeighbours();
	for(int i=0;i<t.getEdgeNumber()*2;i++)
	{
		cout<<" "<<vertexNeighbours[i];
	}
	cout<<endl;

	cout<<"vertexEdgeNumbers:"<<endl;
	swInt *vertexEdgeNumbers=t.getVertexEdgeNumbers();
	for(int i=0;i<t.getVertexNumber();i++)
	{
		cout<<" "<<vertexEdgeNumbers[i];
	}
	cout<<endl;

	cout<<"accuVertexEdgeNumbers:"<<endl;
	swInt *accuVertexEdgeNumbers=t.getAccuVertexEdgeNumbers();
	for(int i=0;i<t.getVertexNumber()+1;i++)
	{
		cout<<" "<<accuVertexEdgeNumbers[i];
	}
	cout<<endl;
}

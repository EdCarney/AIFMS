
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#ifndef PCH_H
#define PCH_H
#endif //PCH_H
#include <string>
#include <limits>
#include <cstring>
#include "mex.h"
#include "simstruc.h"
using namespace std;
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1000000
#define y_width 1000000

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/////////////////////////////////////////////////////////////////////////////////////////////

//#include "heap.h"

double LARGE = std::numeric_limits<double>::infinity();

template <class T> 
class Heap 
{ 
 //  private: 
   public: 
       double* heapCost;     // cost values
       T** heapNode;         // pointers to the corresponding map nodes
                             // (or whatever is stored in the heap)

       int parentOfLast;     // stores the index of the parent of the last node
       int tempInd;          // used to help swap nodes
       double tempCost;      // used to help swap nodes
       T*  tempNode;         // used to help swap nodes

       int capacity;         // the number of things this heap can store without
                             // being increased in size
//   public: 

       int indexOfLast;      // the index of the last node in the heap array

//   private: 

       // sets up the memory for the heap
       // assumes a reserve size of heapCapacity
       void buildHeap(int heapCapacity);

       // compares a node n with its parent, and switches them if the parent's
       // cost is more than the node's cost. Repeats if a switch happens.
       // returns the index that the node ended up at
       int bubbleUp(int n);

       // compares a node n with its children, and switches them if a child's cost
       // is less than the node's cost. Repeats if a switch happens.
       void bubbleDown(int n);

       // cleans up and deletes the memory used by the heap
       void deleteHeap();

       // increases the heap size by a factor of two
       void increaseHeapSize();

//   public: 


     // default constructor
     Heap()
     {
       //printf("building a heap\n");
       buildHeap(10000000);
     }

     // constructor that should be used most of the time
     Heap(int heapCapacity)
     {
       //printf("building a heap with %d elements\n", heapCapacity);
       buildHeap(heapCapacity);
     }


     // destructor
     ~Heap()
     {
       //printf("deleting a heap\n");
       deleteHeap();
     }


     // copy constructor
     //     Heap(const Heap &H)




     // add thisNode to the heap
     void addToHeap(T* thisNode, float keyValue);

     // returns a pointer to the node that is on the top of the heap
     T* topHeap();

     // removes the top valued node from the heap and returns a pointer to it
     T* popHeap();


     // removes the particular node from the heap
     // (even if that node is internal to the heap)
     // and then rebalances the heap, returns a pointer
     // to the node that has been removed
     T* removeNodeFromHeap(T* thisNode);

     // updates the position of the particular node within the heap
     // to reflect the new key value
     // (even if that node is internal to the heap)
     // and then rebalances the heap
     // NOTE: this will insert the node if it is not in the heap already
     void updateNodeInHeap(T* thisNode, float keyValue);


     // prints the heap values on the command line
     void printHeap();

     // returns 1 if heap is good, 0 if bad, also prints a command line message
     // this should only be used for debugging as it will really slow down the code
     int checkHeap();

};

/////////////////////////////////////////////////////////////////////////////////////////////

//#include "heap.cpp"
// sets up the memory for the heap
template <class T> 
void Heap<T>::buildHeap(int heapCapacity)
{

  heapCost = (double*)calloc(heapCapacity, sizeof(double));
  heapNode = (T**)calloc(heapCapacity, sizeof(T*));
  indexOfLast = -1;
  parentOfLast = -1;
  tempCost = LARGE;
  tempNode = NULL;
  capacity = heapCapacity;
  int i;
  for(i = 0; i < heapCapacity; i++)
    heapCost[i] = LARGE;
}

// compares a node n with its parent, and switches them if the parent's
// cost is more than the node's cost. Repeats if a switch happens.
// returns the index that the node ended up at
template <class T>  
int Heap<T>::bubbleUp(int n)
{
  tempInd = (n-1)/2;
  while(n != 0 & heapCost[tempInd] > heapCost[n])
  {
     // swap costs 
     tempCost = heapCost[tempInd]; 
     heapCost[tempInd] = heapCost[n];
     heapCost[n] = tempCost;
     
     // swap graph node pointers
     tempNode = heapNode[tempInd];
     heapNode[tempInd] = heapNode[n];
     heapNode[n] = tempNode;
     
     // update graph node heap index values
     heapNode[tempInd]->heapIndex = tempInd;
     heapNode[n]->heapIndex = n;
     
     // get new node and parent indicies
     n = tempInd;
     tempInd = (n-1)/2;
  } 
  return n;  
}

// compares a node n with its children, and switches them if a child's cost
// is less than the node's cost. Repeats if a switch happens.
template <class T>  
void Heap<T>::bubbleDown(int n)
{
  // find child with smallest value
  if(heapCost[2*n+1] < heapCost[2*n+2])
    tempInd = 2*n+1;
  else
    tempInd = 2*n+2; 

  while(n <= parentOfLast & heapCost[tempInd] < heapCost[n])
  {
     // swap costs 
     tempCost = heapCost[tempInd]; 
     heapCost[tempInd] = heapCost[n];
     heapCost[n] = tempCost;
     
     // swap graph node pointers
     tempNode = heapNode[tempInd];
     heapNode[tempInd] = heapNode[n];
     heapNode[n] = tempNode;  
      
     // update graph node heap index values
     heapNode[tempInd]->heapIndex = tempInd;
     heapNode[n]->heapIndex = n;
     
     // get new node and child indicies
     n = tempInd;
     if(heapCost[2*n+1] < heapCost[2*n+2])
       tempInd = 2*n+1;
     else
       tempInd = 2*n+2;   
  }   
}

// add thisNode to the heap, with the key value
template <class T> 
void Heap<T>::addToHeap(T* thisNode, float keyValue)
{
  if(indexOfLast == capacity-1)
  {
    increaseHeapSize();
  }

  if(thisNode->inHeap == false)
  {
    indexOfLast++;
    parentOfLast = (indexOfLast-1)/2;
    heapNode[indexOfLast] = thisNode;
    heapCost[indexOfLast] = keyValue;
    thisNode->heapIndex = indexOfLast;
    bubbleUp(indexOfLast);     
    thisNode->inHeap = true;
  }
}

// returns a pointer to the node that is on the top of the heap
template <class T> 
T* Heap<T>::topHeap()
{
  if(indexOfLast >= 0)
  {
    return heapNode[0];
  }
  return NULL;
}


// removes the top valued node from the heap and returns a pointer to it
template <class T>  
T* Heap<T>::popHeap()
{
  T* oldTopNode = heapNode[0];
  heapNode[0] = heapNode[indexOfLast];
  heapCost[0] = heapCost[indexOfLast];
  heapNode[0]->heapIndex = 0;
  heapNode[indexOfLast] = NULL;
  heapCost[indexOfLast] = LARGE;
  indexOfLast--;
  parentOfLast = (indexOfLast-1)/2;
  bubbleDown(0);
  oldTopNode->inHeap = false;
  oldTopNode->heapIndex = -1;
  return oldTopNode;
}

// removes the particular node from the heap
// (even if that node is internal to the heap)
// and then rebalances the heap, returns a pointer
// to the node that has been removed
template <class T>  
T* Heap<T>::removeNodeFromHeap(T* thisNode)
{
  if(!thisNode->inHeap)
  {
    return NULL;
  }

  int ind = thisNode->heapIndex;

  // replace the node to be removed from the heap with the last node in the heap
  heapNode[ind] = heapNode[indexOfLast];
  heapCost[ind] = heapCost[indexOfLast];
  heapNode[ind]->heapIndex = ind;

  // set the last position of the heap to have NULL properties
  heapNode[indexOfLast] = NULL;
  heapCost[indexOfLast] = LARGE;

  // decrement the index of the last node
  indexOfLast--;

  // update the parent of the last node
  parentOfLast = (indexOfLast-1)/2;

  // bubble up and down the (formerly) last node in its new position
  ind = bubbleUp(ind);
  bubbleDown(ind);

  thisNode->inHeap = false;
  thisNode->heapIndex = -1;
  return thisNode;
}


// updates the position of the particular node within the heap
// to reflect the new key value
// (even if that node is internal to the heap)
// and then rebalances the heap)
// NOTE: this will insert the node if it is not in the heap already
template <class T>  
void Heap<T>::updateNodeInHeap(T* thisNode, float keyValue)
{
  // we'll just do the easy way
  if(thisNode->inHeap)
  {
    removeNodeFromHeap(thisNode);
  }
  addToHeap(thisNode, keyValue);
}

// prints the heap values on the command line
// note this will break if T does not have field "id"
template <class T> 
void Heap<T>::printHeap()
{
  //printf("heap costs:\n");

  int i,p; 
  char ch[10];
  for(i = 0, p = 0; i <= indexOfLast; i++)
  {    
    //printf("%f ", heapCost[i]);
    if(i == p)
    {
       //printf("\n");
       p = p+2^p;
    }
  }
  //printf("\n\n");
  

  //printf("heap node ids:\n");
  for(i = 0, p = 0; i <= indexOfLast; i++)
  {    
    //printf("(%d) ", heapNode[i]->id);
    if(i == p)
    {
       //printf("\n");
       p = p+2^p;
    }
  }
  //printf("\n\n");
}


// returns 1 if heap is good, 0 if bad, also prints a command line message
template <class T> 
int Heap<T>::checkHeap()
{
  int i;
  for(i = 0; i <= indexOfLast; i++)
  {
    if(heapCost[i] < heapCost[(i-1)/2] || heapNode[i]->heapIndex != i)
    {
      //printf("There is a problem with the heap \n");
      getchar();
      return false;
    }
  } 
  //printf("The heap is OK \n");
  return true;  
}

// cleans up and deletes the memory used by the heap
// note this is called automatically be the destructor
template <class T> 
void Heap<T>::deleteHeap()
{
  free(heapCost);
  heapCost = NULL;
  free(heapNode);
  heapNode = NULL;
  indexOfLast = -1;
  parentOfLast = -1;
  tempInd = -1;
  tempNode = NULL;
  tempCost = LARGE;
}

// increases the heap size by a factor of two
// this should only be used for debugging as it will really slow down the code
template <class T> 
void Heap<T>::increaseHeapSize()
{

  //printf("growing heap\n");
  int oldCapacity = capacity;
  capacity = capacity * 2;
  if(capacity < 100)
  {
    capacity = 100;
  }
  
  double* newHeapCost = (double*)calloc(capacity, sizeof(double));
  T** newHeapNode = (T**)calloc(capacity, sizeof(T*)); 

  std::memcpy(newHeapCost, heapCost, oldCapacity*sizeof(double));
  std::memcpy(newHeapNode, heapNode, oldCapacity*sizeof(T*));

  free(heapCost);
  free(heapNode);

  heapCost = newHeapCost;
  heapNode = newHeapNode;

  //printf("heap can now hold up to %d items\n", capacity);
}

struct Node;        // NOTE: we'll define this in detail later, but we need it now

// make the basic edge data stucture 
struct Edge 
{ 
  Node* startNode;     // the edge starts at this node
  Node* endNode;       // the edge goes to this node

  double edgeCost;  // going along the edge cost this much
}; 
 

// this is a basic graph node for a 2D graph
struct Node
{
  int id;     // this is the id of the node
              // it is alos the position in the node array in the
              // graph data structure where this node is stored

  double x;   // physical x location
  double y;   // physical y location
  double z;   // physical z location

  double costToStart;	// cost from start node to this node
  double heuristic;		// heuristic for A*


  // NOTE, in an undirect graph outgoingEdges and incommingEdges are the same
  // while in a directed graph they are not

  int numOutgoingEdges;  // the number of outgoing edges

  Edge** outgoingEdges;  // these are edges to {neighbors accessible from this node}
                         // outgoingEdges[j] is a pointer to the j-th one
 
  int numIncommingEdges; // the number of outgoing edges

  Edge** incommingEdges; // these are edges from {nodes which this node can be accessed}
                         // incommingEdges[k] is a pointer the k-th one


  Node* parentNode;      // used for graph search, this is a pointer to this node's parent
                         // node in the search tree

  int status;            // used for grah search, 0 = not visited yet
                         //                       1 = in open list (i.e., in heap)
                         //                       2 = in closed list

  // the following fields are assumed by the heap data structure that I've been using

  bool inHeap;           // init to false and set to true iff the node is in the heap
  int heapIndex;         // enables quick access to the node in the heap

};

double calcHeur(Node* thisNode, Node* goalNode)
{
	double newHeuristic;
	newHeuristic = sqrt(pow(thisNode->x - goalNode->x, 2)
		+ pow(thisNode->y - goalNode->y, 2) + pow(thisNode->z - goalNode->z, 2));

	thisNode->heuristic = newHeuristic;
	return newHeuristic;
}

void expand(Heap<Node> &H, Node* thisNode, Node* goalNode)
{
	for (int n = 0; n < thisNode->numOutgoingEdges; n++)
	{
		Edge* thisEdge = thisNode->outgoingEdges[n];  // pointer to this edge
		Node* neighborNode = thisEdge->endNode;  // pointer to the node on the 
												 // other end of this edge
		if (neighborNode->status == 0)  // neighbor has not yet been visited 
		{
			// calculate heuristic score as needed
			// here heuristic is euclidean distance
			neighborNode->heuristic = (neighborNode == goalNode) ? 0 : calcHeur(neighborNode, goalNode);

			// add the neighbor to the heap with key equal to
			// the cost to start plus the node heuristic
			double neighborKey = thisEdge->edgeCost + thisNode->costToStart;	// initially set key to the neighbor node cost to start
			neighborNode->costToStart = neighborKey;							// use this key to set the neighbor node costToStart
			H.addToHeap(neighborNode, neighborKey + neighborNode->heuristic);	// update neighbor node key to costToStart plus heuristic

			// remeber this node as its parent
			neighborNode->parentNode = thisNode;

			// make sure it is in the open list
			neighborNode->status = 1;
		}
		else if (neighborNode->status == 1)  // neighbor is in the open list
		{                                   // (but not the closed list)

		  // will update the parent to this node if the neighbor node has higher
		  // cost to strt otherwise
			if (neighborNode->costToStart > thisEdge->edgeCost + thisNode->costToStart)
			{

				double newNeighborKey = thisEdge->edgeCost + thisNode->costToStart;
				neighborNode->costToStart = newNeighborKey;
				H.updateNodeInHeap(neighborNode, newNeighborKey + neighborNode->heuristic);

				// remeber this node as its parent
				neighborNode->parentNode = thisNode;
			}
		}
	}

	thisNode->status = 2;    // now this node is in the closed list
}

/////////////////////////////////////////////////////////////////////////////////////////////
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void A_Star_S_Func_2_Outputs_wrapper(const real_T *nodeArray,
			const real_T *edgeArray,
			const real_T *numNodesPtr,
			const real_T *numEdgesPtr,
			const real_T *startNodePtr,
			const real_T *goalNodePtr,
			real_T *finalPath,
			real_T *finalPathCost,
			real_T *numNodesCheck,
			real_T *numEdgesCheck)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
class Graph 
	{ 
	public:
		int numNodes;    // the number of nodes in the graph

		Node* nodes;     // an array that contains all of the nodes
						// nodes[i] contains the node with id i

		int numEdges;    // the number of edges in the graph

		Edge* edges;     // an array that contains all of the edges
						// edge[j] contains the j-th edge


		// default constructor
		Graph()
		{
		//printf("building a default graph\n");
		numNodes = 0;
		nodes = NULL;
		numEdges = 0;
		edges = NULL;
		}

		// default destructor
		~Graph()
		{
		//printf("deleting a graph\n");
		if(nodes != NULL)
		{
			for(int n = 0; n < numNodes; n++)
			{
			if(nodes[n].outgoingEdges != NULL)
			{
				free(nodes[n].outgoingEdges);
			}
			nodes[n].numOutgoingEdges = 0;

			if(nodes[n].incommingEdges != NULL)
			{
				free(nodes[n].incommingEdges);
			}
			nodes[n].numIncommingEdges = 0;
			}
			free(nodes);
		}
		numNodes = 0;
		nodes = NULL;

		if(edges != NULL)
		{
			free(edges);
		}
		numEdges = 0;
		edges = NULL;
		}

		bool readGraphFromArrays(const double* nodeArray, const double* edgeArray, int numNodes, int numEdges)
		{          
			nodes = (Node*)calloc(numNodes, sizeof(Node));

			// allocate space for the nodes 
			int id, tempCount = 0;
			float x, y, z;
			for(int n = 0; n < numNodes; n++)
			{
				// allocate the node 
				nodes[n].id = (int) nodeArray[tempCount++] - 1;   // NOTE: switches to c++ 0-based indexing!!!
				nodes[n].x = nodeArray[tempCount++];
				nodes[n].y = nodeArray[tempCount++];
                nodes[n].z = nodeArray[tempCount++];

				nodes[n].numOutgoingEdges = 0;  // dummy value for now, reset lower down
				nodes[n].outgoingEdges = NULL;  // dummy value for now, reset lower down
				nodes[n].numIncommingEdges = 0; // dummy value for now, reset lower down
				nodes[n].incommingEdges = NULL; // dummy value for now, reset lower down

				nodes[n].parentNode = NULL;     // used for graph search
				nodes[n].status = 0;            // everything starts not visited

				nodes[n].inHeap = false;        // used by heap
				nodes[n].heapIndex = -1;        // used by heap
			}
            
			// allocate space for the edges
			edges = (Edge*)calloc(numEdges, sizeof(Edge));

			// now read the edges one at a time
			int startNodeID, endNodeID;
			float cost;
			tempCount = 0;
			for(int m = 0; m < numEdges; m++)
			{
				startNodeID = (int) edgeArray[tempCount++];
				endNodeID = (int) edgeArray[tempCount++];
				cost = edgeArray[tempCount++];

				// allocate the edge 
				// NOTE: switches to c++ 0-based indexing!!!
				// also note that edges store pointers to nodes that they are the edge for
				edges[m].startNode = &nodes[startNodeID-1];
				edges[m].endNode = &nodes[endNodeID-1];
                
				// cost from file
				edges[m].edgeCost = cost;
                
			}

			// now we need to go through the nodes and set up the neighbor
			// lists for each node. In this code we are using edges arrays
			// for this part. we will use three passes, on the first
			// pass we figure out how many incomming and outgoing edges
			// each node has, next we allocte all required space,
			// finally, we actually set the edges up


			// we'll also init some temp variables to help us
			// remember how many edges have been set so far in the final pass
			int* outgoingEdgesCount = (int*)calloc(numNodes, sizeof(int));
			int* incommingEdgeCount = (int*)calloc(numNodes, sizeof(int));
			for(int n = 0; n < numNodes; n++)
			{
				outgoingEdgesCount[n] = 0;
				incommingEdgeCount[n] = 0; 
				nodes[n].numOutgoingEdges = 0;   // set this to 0 while we're at it
				nodes[n].numIncommingEdges = 0;  // set this to 0 while we're at it
			}   
            
			// first pass: go through all edges to see how many incomming 
			// and outgoing edges each node has, and recored these values      
			for(int m = 0; m < numEdges; m++)
			{
				outgoingEdgesCount[edges[m].startNode->id]++;
				incommingEdgeCount[edges[m].endNode->id]++;
			}      
            
			// second: go through all nodes and allocate space for their
			// neighbor lists, 
			for(int n = 0; n < numNodes; n++)
			{
				nodes[n].outgoingEdges = (Edge**)calloc(outgoingEdgesCount[n], sizeof(Edge*));
				nodes[n].incommingEdges = (Edge**)calloc(incommingEdgeCount[n], sizeof(Edge*));
			}
            
			// third: go through all the edges again and record pointers to them 
			// in the neighbor lists
			for(int m = 0; m < numEdges; m++)
			{
				// use local pointers to the start and end nodes to make 
				// things things easier to read
				Node* startNodePtr = edges[m].startNode; 
				Node* endNodePtr = edges[m].endNode; 
				//printf("start node: %d\n", startNodePtr->id);
				//printf("end node: %d\n", endNodePtr->id);

				startNodePtr->outgoingEdges[startNodePtr->numOutgoingEdges] = &edges[m];
				endNodePtr->incommingEdges[endNodePtr->numIncommingEdges] = &edges[m];

				startNodePtr->numOutgoingEdges++;
				endNodePtr->numIncommingEdges++;
			}
            
			// cleanup
			free(outgoingEdgesCount);
			free(incommingEdgeCount);
            
            return true;
        }
	};
    
//int numNodes = *((int*)numNodesPtr);
//int numEdges = *((int*)numEdgesPtr);
  
int numNodes = (int)numNodesPtr[0];
int numEdges = (int)numEdgesPtr[0];
int startNodeIndex = (int)startNodePtr[0] - 1;  // update indexing!
int goalNodeIndex = (int)goalNodePtr[0] - 1;    // update indexing!
    
Graph G;
G.readGraphFromArrays(nodeArray, edgeArray, numNodes, numEdges);
Heap<Node> H(numNodes);

// these are pointers to the start and end nodes
Node* startNode = &G.nodes[startNodeIndex];
Node* goalNode = &G.nodes[goalNodeIndex];

// initialize the start node
// the key is equal to the costToStart (zero for start) plus the heuristic score

startNode->status = 1;									// now the start node is in the open list
startNode->costToStart = 0;								// set costToStart
startNode->heuristic = calcHeur(startNode, goalNode);	// set heuristic
double key = startNode->heuristic;						// use heuristic to set key (costToStart zero for startNode)
H.addToHeap(startNode, key);

// do while there are nodes left in the heap
// this will stop once the goal node is reached,
// but will still do the last expansion with the goal node
while (H.topHeap() != goalNode && H.topHeap() != NULL)	// add a check if the heap is empty
	{
        Node* thisNode = H.popHeap();
		    expand(H, thisNode, goalNode);
	}

if (H.topHeap() != NULL)								// if heap empty then no path possible; print warning message
{
    Node* thisNode = goalNode;
    int iterator = 0;
    while(thisNode != NULL)
    {
        // format is id, x, y    
        // NOTE: incrimenting ids by 1 to convert to 1-based-indexing 
        // fprintf(pFile, "%d, %f, %f\n",thisNode->id+1, thisNode->x, thisNode->y);
        finalPath[iterator++] = thisNode->id+1;
        thisNode = thisNode->parentNode;
    }
    finalPathCost[0] = (float) goalNode->costToStart;
}

numNodesCheck[0] = numNodes;
numEdgesCheck[0] = numEdges;
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}



//                              -*- Mode: C++ -*- 
// Graph.cc
// Copyright © 2002 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Mon Feb 18 16:07:09 2002
// Last Modified By : 
// Last Modified On : 
// Update Count     : 0
// Status           : Unknown.
// 



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream.h>
#include <iomanip.h>
#include <values.h>

#include "Graph.h"


template< class Node, class Edge >
Graph< Node, Edge >::Graph (const Graph &other)
  : adjgraph (other)
{
  if (other.mNode.size() > 0) 
    {
      mNodes = other.mNodes;  
      mNodeValues = other.mNodeValues;
      mEdges = other.mEdges;  
      mEdgeValues = other.mEdgeValues;
      mEdgeNodes = other.mEdgeNodes;
    }
}


template< class Node, class Edge >
Graph< Node, Edge >::~Graph ()
{
  mNodes.erase (mNodes.begin(), mNodes.end());
  mNodeValues.erase (mNodeValues.begin(), mNodeValues.end());
  mEdges.erase (mEdges.begin(), mEdges.end());
  mEdgeValues.erase (mEdgeValues.begin(), mEdgeValues.end());
  mEdgeNodes.erase (mEdgeNodes.begin(), mEdgeNodes.end());
  erase (begin(), end()); 
}


template< class Node, class Edge >
Graph< Node, Edge >& Graph< Node, Edge >::operator= (const Graph &other)
{
  if  (this != &other) { 
    if  (other.mNodes.size() > 0) 
      {
	adjgraph::operator= (other);
	mNodes = other.mNodes;  
	mNodeValues = other.mNodeValues;
	mEdges = other.mEdges;  
	mEdgeValues = other.mEdgeValues;
	mEdgeNodes = other.mEdgeNodes;
      }
  }
  return *this;
}


template< class Node, class Edge >
const Node &Graph< Node, Edge >::getNode (int index) const
{
  assert (index >= 0 && index < (int)mNodes.size());
  return mNodes[ index ]; 
}


template< class Node, class Edge >
Node &Graph< Node, Edge >::getNode (int index)
{
  assert (index >= 0 && index < (int)mNodes.size());
  return mNodes[ index ]; 
}

  
template< class Node, class Edge >
int Graph< Node, Edge >::getNodeValue (int index) const
{
  if  (index < 0 || index > (int)mNodeValues.size()) {
    return MAXINT;
  }
  return mNodeValues[ index ]; 
}



template< class Node, class Edge >
void Graph< Node, Edge >::setNodeValue (int v)
{
  assert (index >= 0 && index < (int)mNodeValues.size());
  mNodeValues[ index ] = v;
}

 
template< class Node, class Edge >
int Graph< Node, Edge >::getNodeIndex (const Node &n) const
{
  vector< Node >::const_iterator i = ::find (mNodes.begin(), mNodes.end(), n);
  if  (i == mNodes.end()) return -1;
  return  (i - mNodes.begin());
}
 
 
template< class Node, class Edge >
const Edge &Graph< Node, Edge >::getEdge (int index) const
{
  assert (index >= 0 && index < (int)mEdge.size());
  return mEdges[ index ]; 
}


template< class Node, class Edge >
Edge &Graph< Node, Edge >::getEdge (int index)
{
  assert (index >= 0 && index < (int)mEdge.size());
  return mEdges[ index ]; 
}


template< class Node, class Edge >
const Edge &Graph< Node, Edge >::getEdge (int a, int b) const
{
  adjgraph::const_iterator i = (*this).find (a);
  if  (i == end()) {
    cerr << "Index out of range" << endl;
    exit (EXIT_FAILURE);
  }
  adjlist::const_iterator j = i->second.find (b);
  if  (j == i->second.end()) {
    cerr << "Index out of range" << endl;
    exit (EXIT_FAILURE);
  }
  return mEdges[ j->second ];
}


template< class Node, class Edge >
Edge &Graph< Node, Edge >::getEdge (int a, int b)
{
  adjgraph::const_iterator i = (*this).find (a);
  if  (i == end()) {
    cerr << "Index out of range" << endl;
    exit (EXIT_FAILURE);
  }
  adjlist::const_iterator j = i->second.find (b);
  if  (j == i->second.end()) {
    cerr << "Index out of range" << endl;
    exit (EXIT_FAILURE);
  }
  return mEdges[ j->second ];
}


template< class Node, class Edge >
const pair< int, int > Graph< Node, Edge >::getNodes (const Edge &e) const
{
  vector< Edge >::const_iterator i = ::find (mEdges.begin(), mEdges.end(), e);
  if  (i == mEdges.end()) return make_pair (-1, -1);
  return mEdgeNodes[(i-mEdges.begin ())];
}

template< class Node, class Edge >
int Graph< Node, Edge >::getEdgeValue (int index) const
{
  if  (index < 0 || index > (int)mEdges.size()) {
    return MAXINT;
  }
  return mEdgeValues[ index ]; 
}


template< class Node, class Edge >
void Graph< Node, Edge >::setEdgeValue (int index, int v)
{
  assert (index >= 0 && index < (int)mEdgeValues.size ());
  mEdgeValues[ index ] = v; 
}


template< class Node, class Edge >
int Graph< Node, Edge >::getEdgeIndex (const Edge &e) const
{
  vector< Edge >::const_iterator i = ::find (mEdges.begin(), mEdges.end(), e);
  if  (i == mEdge.end()) return -1;
  return  (i - mEdge.begin());
}


template< class Node, class Edge >
int Graph< Node, Edge >::getEdgeIndex (int a, int b) const
{
  adjgraph::const_iterator i = (*this).find (a);
  if  (i == end()) return -1;
  adjlist::const_iterator j = i->second.find (b);
  if  (j == i->second.end()) return -1;
  return j->second;
}


// Check if there is a directed edge between a and b
// Returns false even if there is a directed edge between b and a
template< class Node, class Edge >
bool Graph< Node, Edge >::isConnected (int a, int b) const // a=origin, b=destination
{ 
  adjgraph::const_iterator i = (*this).find (a);
  if  (i == end()) return false;
  adjlist::const_iterator j = i->second.find (b);
  if  (j == i->second.end()) return false;
  return true;
}

  
template< class Node, class Edge >
int Graph< Node, Edge >::addNode (const Node &n, const int v = 0)
{
  // TODO : check redundancy
  unsigned int size = mNodes.size();
  mNodes.push_back (n);
  mNodeValues.push_back (v);
  return size;
}

template< class Node, class Edge >
int Graph< Node, Edge >::addEdgeById (int a, int b, const Edge &e, 
				      bool oriented = true, const int v = 0)
{
  // TODO : Careful about redundancy...
  unsigned int size = mEdges.size();
  
  mEdges.push_back (e);
  mEdgeNodes.push_back (make_pair (a,b));
  mEdgeValues.push_back (v);
  
  (*this)[ a ][ b ] = size;
  if  (!oriented) {
    (*this)[ b ][ a ] = size;
  }
  return size;
}
  
template< class Node, class Edge >
int Graph< Node, Edge >::addEdge (const Node &a, const Node &b, const Edge &e, 
				  bool oriented = true, const int v = 0)
{
  int ida = getNodeIndex (a);
  int idb = getNodeIndex (b);
  return addEdgeById (ida, idb, e, oriented, v);
}
    
template< class Node, class Edge >
vector< Path< int > > Graph< Node, Edge >::shortest_path (int root) const
{
  int i, j;
  int n = (int)size();
  vector< Path< int > > P (n);    // path description
  vector< int > C;                  // node set
  vector< int >::iterator k;
  
  // Initialize ---
  for  (i=0; i<n; i++) {
    if  (i == root) { 
      P[ i ].setValue (0);
    } else {
      C.push_back (i);
      if  (isConnected (root, i)) {
	P[ i ].setValue (getEdgeValue (getEdgeIndex (root, i)));
	P[ i ].push_back (root); 
	P[ i ].push_back (i);
      } else {
	P[ i ].setValue (MAXINT);
      }
    }
  }
  
  //        cout << "P =  ("; for  (i=0; i<n; i++) { cout << P[ i ] << " "; } cout << ")" << endl;
  //        cout << "C =  ("; for  (i=0; i<C.size(); i++) { cout << C[ i ] << " "; } cout << ")" << endl;
  
  // Execute ---
  for  (i=0; i<n-2; i++) {
    vector< int >::iterator min_iter = C.begin();
    int min_value = P[ *min_iter ].getValue(); // in C
    int min_index;
    
    for  (k=C.begin(); k!=C.end(); k++) {
      if  (P[ *k ].getValue() < min_value) {
	min_value = P[ *k ].getValue();
	min_iter = k;
      } 
    }
    min_index = *min_iter;
    C.erase (min_iter);
    
    if  (min_value == MAXINT) break; // in case there is no better element...
    
    for (k=C.begin(); k!=C.end(); k++) {
      int old_val = P[ *k ].getValue();
      int new_val;
      if (isConnected (min_index, *k))
	new_val = min_value + getEdgeValue (getEdgeIndex (min_index, *k));
      else 
	new_val = MAXINT;
      if  (old_val > new_val) {
	P[ *k ] = P[ min_index ];
	P[ *k ].push_back (*k);
	P[ *k ].setValue (new_val);
      }
    }
  }
  return P;
}
  
template< class Node, class Edge >
vector< Edge > Graph< Node, Edge >::minimum_spanning_tree (void) const
{
  vector< int > nearest (size (), -1);
  vector< int > mindist (size (), MAXINT);
  vector< Edge > T;
  int i, j;

  for (i=1; i<(int)size (); ++i) {
    nearest[i] = 0;
    mindist[i] = getEdgeValue (getEdgeIndex (i, 0));
  }

  for (int x=0; x<(int)size ()-1; ++x) {
    int min = MAXINT;
    int k = -1;
    for (j=1; j<(int)size (); ++j) {
      if (mindist[j] >= 0 && mindist[j]<min) {
	min = mindist[j];
	k = j;
      }
    }

    // This is a test to see if we stay in the same connected component
    if (k!=-1) 
      {  
	T.push_back (getEdge (nearest[k], k));
	
	mindist[k] = -1;
	
	for (j=1; j<(int)size (); ++j) {
	  int val = getEdgeValue (getEdgeIndex (j, k));
	  if (val < mindist[j]) {
	    mindist[j] = val;
	    nearest[j] = k;
	  }
	}
      }
  }
  
  return T;
}

// I/O  -----------------------------------------------------------------


template< class Node, class Edge >
ostream &Graph< Node, Edge >::output (ostream &out) const 
{
  if  (size() == 0) return out;
  unsigned int i, j;
  int id;
  
  for  (i=0; i<mNodes.size(); i++) {
    out << "Node[ " << i << " ]: " << mNodes[ i ] 
	<< " (value: " << mNodeValues[ i ] << ") " << endl;
  }
  out << endl;
  for  (j=0; j<mEdges.size(); j++) {
    out << "Edge[ " << j << " ]: " << mEdges[ j ] 
	<< " (value: " << mEdgeValues[ j ] << ") " << endl; 
  }
  out << endl;
  
  out << "Matrix:" << endl;
  for (i=0 ; i<mNodes.size() ; i++) {
    for (j=0 ; j<mNodes.size() ; j++) {
      id = getEdgeIndex (i, j);
      if  (id != -1)
	out << setw (2) << id << " ";
      else
	out << " - ";
    }
    out << "\n";
  }
  return out;
}

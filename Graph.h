//                              -*- Mode: C++ -*- 
// Graph.h
// Copyright © 2002 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Mon Feb 18 16:07:09 2002
// Last Modified By : 
// Last Modified On : 
// Update Count     : 0
// Status           : Unknown.
// 


#ifndef _Graph_h_
#define _Graph_h_

#include <iostream.h>
#include <map.h>

/**
 * @short Description
 *
 * Long Description
 *
 * @author Patrick Gendron
 */
template< class Node >
class Path : public vector< Node >
{
  int mValue;
  
public:
  // LIFECYCLE ------------------------------------------------------------

  /**
   * Initializes the object.
   */
  Path () : vector< Node >() { mValue = 0; }

  /**
   * Initializes the object with the other's content.
   * @param other the object to copy.
   */
  Path (const Path &other) : vector< Node > (other) { mValue = other.mValue; }

  /**
   * Destroys the object.
   */
  ~Path () {}

  // OPERATORS ------------------------------------------------------------

  /**
   * Assigns the object with the other's content.
   * @param other the object to copy.
   * @return itself.
   */
  Path& operator= (const Path &other)
  {
    if  (this != &path) {
      vector< Node >::operator= (path);
      mValue = path.mValue;
    } 
    return *this;
  }

  bool operator< (const Path &path) const { // used in sort
    return  (size () < path.size ());
  }

  void setValue (int value) { mValue = value; }
  int getValue (void) const { return mValue; }

  // ACCESS ---------------------------------------------------------------

  // METHODS --------------------------------------------------------------

  // I/O  -----------------------------------------------------------------

};

template< class Node >
ostream &operator<< (ostream &out, const Path< Node > &path) 
{
  out << "[ " << flush;
  for  (unsigned int i=0; i<path.size(); i++) 
    out << path[ i ] << " ";
  out << "] = " << path.mValue << flush;
  return out;
}



/**
 * @short Description
 *
 * Long Description
 *
 * @author Patrick Gendron
 */
template< class Node, class Edge >
class Graph : public map< int, map< int, int > >
{
public:

typedef map< int, int > adjlist;    // Adjacency List
typedef map< int, adjlist > adjgraph;  // Adjacency Graph  

protected:

  vector< Node > mNodes;
  vector< int > mNodeValues;
  
  vector< Edge > mEdges;
  vector< int >  mEdgeValues;
  vector< pair< int, int > > mEdgeNodes;

public:

  // LIFECYCLE ------------------------------------------------------------

  /**
   * Initializes the object.
   */
  Graph () {}

  /**
   * Initializes the object with the other's content.
   * @param other the object to copy.
   */
  Graph (const Graph &other);

  /**
   * Destroys the object.
   */
  ~Graph ();

  // OPERATORS ------------------------------------------------------------

  /**
   * Assigns the object with the other's content.
   * @param other the object to copy.
   * @return itself.
   */
  Graph& operator= (const Graph &other);

  // ACCESS ---------------------------------------------------------------

  size_t size () const { return mNodes.size (); }

  // METHODS --------------------------------------------------------------

  const vector< Node > &getNodes (void) const { return mNodes; }
  vector< Node > &getNodes (void) { return mNodes; }
  
  const Node &getNode (int index) const;
  Node &getNode (int index);
  
  int getNodeValue (int index) const;
  void setNodeValue (int v);
  
  int getNodeIndex (const Node &n) const;
  
  const vector< Edge > &getEdges (void) const { return mEdges; }
  vector< Edge > &getEdges (void) { return mEdges; }

  const vector< pair< int, int > > &getEdgeNodes (void) const 
  { return mEdgeNodes; }

  const Edge &getEdge (int index) const;
  Edge &getEdge (int index);

  const Edge &getEdge (int a, int b) const;
  Edge &getEdge (int a, int b); // a=origin, b=destination

  const pair< int, int > getNodes (const Edge &e) const;

  int getEdgeValue (int index) const;
  void setEdgeValue (int index, int v);

  int getEdgeIndex (const Edge &e) const; 
  int getEdgeIndex (int a, int b) const; // a=origin, b=destination

  // Check if there is a directed edge between a and b
  // Returns false even if there is a directed edge between b and a
  bool isConnected (int a, int b) const; // a=origin, b=destination
  
  int addNode (const Node &n, const int v = 0);
  int addEdgeById (int a, int b, const Edge &e, 
		   bool oriented = true, const int v = 0);
  
  int addEdge (const Node &a, const Node &b, const Edge &e, 
	       bool oriented = true, const int v = 0);

    
  vector< Path< int > > shortest_path (int root) const;
  
  vector< Edge > minimum_spanning_tree (void) const;

  // I/O  -----------------------------------------------------------------
  ostream &output (ostream &out) const;

};

template< class Node, class Edge >
ostream &operator<< (ostream &out, const Graph< Node, Edge > &gr)
{
  return gr.output (out);
}

#include "Graph.cc"

#endif

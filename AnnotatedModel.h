//                              -*- Mode: C++ -*- 
// AnnotatedModel.h
// Copyright © 2001, 2002 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Fri Nov 16 13:46:22 2001
// Last Modified By : 
// Last Modified On : 
// Update Count     : 0
// Status           : Unknown.
// 


#ifndef _AnnotatedModel_h_
#define _AnnotatedModel_h_

#include <vector.h>
#include <map.h>

#include "mccore/CResidue.h"
#include "mccore/CResId.h"
#include "mccore/Model.h"
#include "mccore/ResidueType.h"
#include "mccore/McCore.h"

#include "mcpl/mcpl.h"
#include "mcpl/Relation.h"
#include "mcpl/Conformation.h"
#include "mcpl/PairingPattern.h"
#include "mcpl/PropertyType.h"
#include "mcpl/MaximumFlow.h"
#include "mcpl/PairAnnote.h"

#include "Graph.h"

typedef int node;
typedef int edge;
typedef int strandid;

enum stype { BULGE_OUT, BULGE, INTERNAL_LOOP, LOOP, OTHER };


struct Helix : public vector< pair< node, node > > {
};

struct Strand : public pair< node, node > {
  stype type;
  strandid ref;
};

/**
 * @short Description
 *
 * Long Description
 *
 * @author Patrick Gendron
 */
class AnnotatedModel
{
  Model *model;

  /**
   * Conformations:  in sequences order.
   */
  vector< Conformation > conformations;

  vector< int > sequence_mask;
  vector< int > sequence_length;

  vector< int > helix_mask;
  vector< Helix > helices;

  vector< int > strand_mask;
  vector< Strand > strands;

  vector< int > tertiary_mask;

  /**
   * Relations
   */
  vector< Relation > relations;

  Graph< node, edge > graph;
  vector< char > marks;

  map< CResId, CResId > translation;

  int nb_pairings;
  int nb_connect;

public:

  // LIFECYCLE ------------------------------------------------------------

  /**
   * Initializes the object.
   */
  AnnotatedModel (Model *m);

  /**
   * Destroys the object.
   */
  ~AnnotatedModel ();

  // OPERATORS ------------------------------------------------------------

  // ACCESS ---------------------------------------------------------------

//    vector< Sequence > & getSequences () { return sequences; }
//    Sequence & getSequence (int i) { return sequences[i]; }
  
  CResId & getResId (node i) { 
    return translation[(const CResId&)(conformations[i].getRes ())]; 
  }
  CResId & getRefId (node i, node j) { 
    return translation[(const CResId&)(relations[graph.getEdge(i, j)].getRef ())]; 
  }
  CResId & getResId (node i, node j) { 
    return translation[(const CResId&)(relations[graph.getEdge(i, j)].getRes ())]; 
  }
  
  t_Residue* getType (node i) {
    static t_Residue* unknown = NULL;
    //if (!unknown) unknown = new rt_Misc ("X");
    t_Residue* t = conformations[i].getRes ().GetType ();
    if (t->is_A ()) return r_A;
    else if (t->is_C ()) return r_C;
    else if (t->is_G ()) return r_G;
    else if (t->is_T ()) return r_T;
    else if (t->is_U ()) return r_U;
    else return unknown;
  }

  bool isPairing (edge e) {
    return relations[e].isPairing ();
  }
  
  bool isHelixPairing (edge e) {
    return (relations[e].is (p_XX) || relations[e].is (p_XIX) ||
	    relations[e].is (p_XXVIII));
  }

  bool isAdjacent (edge e) {
    return relations[e].isAdjacent ();
  }

  bool isStacking (edge e) {
    return relations[e].isStacking ();
  }

  bool isPairing (node i, node j) {
    if (graph.isConnected (i, j))
      return relations[graph.getEdge (i, j)].isPairing ();
    return false;
//      adjgraph::iterator gi;
//      adjlist::iterator gj;
//      gi = graph.find (i);
//      if (gi == graph.end ()) return false;
//      gj = gi->second.find (j);
//      if (gj == gi->second.end ()) return false;
//      return relations[gj->second].isPairing ();
  }

  // METHODS --------------------------------------------------------------

  /**
   * Creates the annotation and determines the sequences present in the 
   * structure.
   */
  void annotate ();

  /**
   * Finds all helices.
   */
  void findHelices ();

  /**
   * 
   */
  void findBulges ();

  /**
   *
   */
  void findStrands ();

  /**
   *
   */
  void findKissingHairpins ();

  /**
   * 
   */
  void findPseudoknots ();

  // I/O  -----------------------------------------------------------------

  void dumpSequences ();
  void dumpGraph ();
  void dumpConformations ();
  void dumpPairs ();
  
  void dumpHelices ();
  void dumpTriples ();
  void dumpStrands ();
  void dumpStacks ();
  void classifyStrands ();
  
  void dumpMcc (const char* pdbname);

};

#endif

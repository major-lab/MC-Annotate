//                              -*- Mode: C++ -*- 
// SStructure.cc
// Copyright © 2006 Laboratoire de Biologie Informatique et Théorique
//                  Université de Montréal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Tue Nov 21 16:24:43 2006
// $Revision$
// $Id$
// 


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "SStructure.h"



using namespace annotate
{

  void
  SStructure::buildSequences ()
  {
    set< AnnotateModel::label > vertexSet;
    AnnotateModel::iterator it;

    for (it = begin (); end () != it; ++it)
      {
	vertexSet.insert (getVertexLabel (&*it));
      }
    while (! vertexSet.empty ())
      {
	AnnotateModel::label loc;
	Sequence sequence;

	vertexSet.erase (loc = *vertexSet.begin ());
	sequence.push_back (loc);
	buildSequence5p (vertexSet, sequence, loc);
	buildSequence3p (vertexSet, sequence, loc);
	if (1 < sequence.size ())
	  {
	    sequences.push_back (sequence);
	  }
      }
  }


  void
  SStructure::buildSequence5p (set< AnnotateModel::label > &vertexSet, Sequence &seq, AnnotateModel::label loc)
  {
    list< AnnotateModel::label > neighbors;
    list< AnnotateModel::label >::const_iterator nit;

    neighbors = internalOutNeighborhood (loc);
    for (nit = neighbors.begin (); neighbors.end () != nit;)
      {
	AnnotateModel::label nloc;

	if (internalGetEdge (loc, nloc = *nit)->is (PropertyType::pAdjacent5p))
	  {
	    vertexSet.erase (nloc);
	    seq.push_back (nloc);
	    loc = nloc;
	    neighbors = internalOutNeighborhood (loc);
	    nit = neighbors.begin ();
	    continue;
	  }
	++nit;
      }
  }
  

  void
  SStructure::buildSequence3p (set< AnnotateModel::label > &vertexSet, Sequence &seq, AnnotateModel::label loc)
  {
    list< AnnotateModel::label > neighbors;
    list< AnnotateModel::label >::const_iterator nit;

    neighbors = internalOutNeighborhood (loc);
    for (nit = neighbors.begin (); neighbors.end () != nit;)
      {
	AnnotateModel::label nloc;

	if (internalGetEdge (loc, nloc = *nit)->is (PropertyType::pAdjacent3p))
	  {
	    vertexSet.erase (nloc);
	    seq.insert (seq.begin (), nloc);
	    loc = nloc;
	    neighbors = internalOutNeighborhood (loc);
	    nit = neighbors.begin ();
	    continue;
	  }
	++nit;
      }
  }


  void
  SStructure::annotate ()
  {
    
  }

}

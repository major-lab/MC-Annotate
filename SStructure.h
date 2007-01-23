//                              -*- Mode: C++ -*- 
// SStructure.h
// Copyright © 2006 Laboratoire de Biologie Informatique et Théorique
//                  Université de Montréal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Tue Nov 21 16:07:20 2006
// $Revision$
// $Id$
// 

#ifndef _annotate_SStructure_h_
#define _annotate_SStructure_h_

#include <iostream>
#include <vector>

#include "mccore/AnnotateModel.h"
#include "mccore/UndirectedGraph.h"

using namespace mccore;
using namespace std;



namespace annotate
{
  
  /**
   * 
   * @author Martin Larose (<a href="larosem@iro.umontreal.ca">larosem@iro.umontreal.ca</a>)
   * @version $Id$
   */
  class SStructure : public UndirectedGraph< AnnotateModel::label, char, char, unsigned int >
  {
    typedef UndirectedGraph< GraphModel::label, char, char, unsigned int > super;

    const AnnotateModel *graph;
    
  public:

    // LIFECYCLE ------------------------------------------------------------

    SStructure (const AnnotateModel &am) : super (), graph (&am) { }

    virtual ~SStructure () { }

    // OPERATORS ------------------------------------------------------------
 
    // ACCESS ---------------------------------------------------------------

    // METHODS --------------------------------------------------------------

    /**
     * Builds the secondary structure on top of the sequence.
     */
    void annotate ();
    
    // I/O  -----------------------------------------------------------------

 };
  
}

#endif

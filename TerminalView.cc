//                              -*- Mode: C++ -*- 
// TerminalView.cc
// Copyright © 2004 Laboratoire de Biologie Informatique et Théorique
//                  Université de Montréal.
// Author           : Martin Larose <larosem@iro.umontreal.ca>
// Created On       : Mon Dec  6 18:11:47 2004
// $Revision$
// $Id$


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "mccore/Messagestream.h"

#include "AnnotateModule.h"
#include "Module.h"
#include "TerminalView.h"

using namespace std;



namespace annotate
{
  void
  TerminalView::show () const
  {
    *ms (0) << "Residue conformations -------------------------------------------" << endl;
    ((AnnotateModule*) module)->dumpConformations ();
    *ms (0) << endl;
    
    ((AnnotateModule*) module)->dumpStacks ();
    *ms (0) << endl;
    
    *ms (0) << "Base-pairs ------------------------------------------------------" << endl;
    ((AnnotateModule*) module)->findKissingHairpins ();
    *ms (0) << endl;
    
    *ms (0) << "Triples ---------------------------------------------------------" << endl;
    ((AnnotateModule*) module)->dumpTriples ();
    *ms (0) << endl;
    
    *ms (0) << "Helices ---------------------------------------------------------" << endl;
    ((AnnotateModule*) module)->dumpHelices ();
    *ms (0) << endl;
    
    *ms (0) << "Strands ---------------------------------------------------------" << endl;
    ((AnnotateModule*) module)->dumpStrands ();
    *ms (0) << endl;
    
    *ms (0) << "Various features ------------------------------------------------" << endl;
    
    *ms (3) << "Finding pseudoknots..." << endl;
    ((AnnotateModule*) module)->findPseudoknots ();
    *ms (0) << endl;
    
    *ms (0) << "Sequences -------------------------------------------------------" << endl;
    ((AnnotateModule*) module)->dumpSequences ();
    *ms (0) << endl;
  }
}

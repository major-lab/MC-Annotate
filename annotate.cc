//                              -*- Mode: C++ -*- 
// annote.cc
// Copyright © 2001, 2002 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Fri May 18 09:38:07 2001
// Last Modified By : 
// Last Modified On : 
// Update Count     : 0
// Status           : Unknown.
// 

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream.h>
#include <stdlib.h>
#include <malloc.h>

#include "mccore/McCore.h"
#include "mccore/zfPdbstream.h"
#include "mccore/CResidue.h"
#include "mccore/Model.h"
#include "mccore/Algo.h"
#include "mccore/CPoint3D.h"
#include "mccore/Messagestream.h"
#include "mccore/zfstream.h"


#include "mcpl/mcpl.h"
#include "mcpl/Relation.h"
#include "mcpl/Conformation.h"
#include "mcpl/PairingPattern.h"
#include "mcpl/PropertyType.h"
#include "mcpl/MaximumFlow.h"
#include "mcpl/PairAnnote.h"

#include "AnnotatedModel.h"


#define NONE_MARK  0
#define HELIX_MARK 1
#define PAIR_MARK  2
#define BULGE_MARK 3
#define LOOP_MARK  4


int main (int argc, char* argv[])
{
  McCoreInit ();
  mcplInit ();

  gOut.setVerboseLevel (1);

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " pdbfile" << endl;
    exit (EXIT_FAILURE);
  }

  izfPdbstream fin (argv[1]);
  Model *model = new Model;

  cerr << "Loading structure... ";

  if (!fin) {
    cerr << "File could not be opened" << endl;
    exit (EXIT_FAILURE);
  }

  fin >> *model;

  Model prot = *model;

  model->keepNucleicAcid ();

  cerr << "found " << model->size () << " nucleic acid residues." << endl;

  prot.keepAminoAcid ();

  cerr << "found " << prot.size () << " amino acid residues." << endl;

  if (model->size () > 0) {

    if (strncmp (argv[0]+strlen(argv[0])-7, "pdb2mcc", 7) == 0) {
      AnnotatedModel amodel (model);
      amodel.annotate ();
      amodel.dumpMcc (argv[1]);
    } else {
      
      //      clock_t startTime, endTime;
      //      startTime = clock();
      
      
      AnnotatedModel amodel (model);
      
      cout << "Annotating..." << flush;
      amodel.annotate ();
      cout << "done." << endl;
      cout << "Finding helices..." << flush;
      amodel.findHelices ();
      cout << "done." << endl;
      cout << "Finding strands..." << flush;
      amodel.findStrands ();
      amodel.classifyStrands ();
      cout << "done." << endl << endl;
      
      //      amodel.dumpGraph ();
      
      amodel.dumpConformations ();
      
      cout << endl;
      
      amodel.dumpPairs ();
      
      cout << endl;
      
      amodel.dumpTriples ();
      
      cout << endl;
      
      amodel.dumpStacks ();
      cout << endl;
      amodel.dumpHelices ();
      cout << endl;
      amodel.dumpStrands ();
      cout << endl;
      amodel.findKissingHairpins ();
      cout << endl;
      amodel.findPseudoknots ();
      cout << endl;
      amodel.dumpSequences ();
      
      cout << "-----------------------------------------------------------------" << endl;
    //      endTime = clock();
    //      cout << "Elapsed time: " << ( (float)(endTime - startTime) / (float)CLOCKS_PER_SEC ) << " sec." << endl;
    }
  }

  return EXIT_SUCCESS;
}






//                              -*- Mode: C++ -*- 
// annote.cc
// Copyright © 2001, 2002, 2003, 2004 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Fri May 18 09:38:07 2001
// Last Modified By : Patrick Gendron
// Last Modified On : Fri Jan  9 10:24:58 2004
// Update Count     : 112
// Status           : Unknown.
// 

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream>
#include <stdlib.h>
#include <malloc.h>
#include <getopt.h>

#include "mccore/zfPdbstream.h"
#include "mccore/fPdbstream.h"
#include "mccore/zfBinstream.h"
#include "mccore/Residue.h"
#include "mccore/Model.h"
#include "mccore/Algo.h"
#include "mccore/Vector3D.h"
#include "mccore/Messagestream.h"
#include "mccore/zfstream.h"
#include "mccore/Relation.h"
#include "mccore/PropertyType.h"
#include "mccore/stlio.h"

#include "AnnotatedModel.h"

#include "MolecularComplex.h"

#define NONE_MARK  0
#define HELIX_MARK 1
#define PAIR_MARK  2
#define BULGE_MARK 3
#define LOOP_MARK  4

int __current_option;
const char* __shortopts = "ae:f:ms:vV";
bool __sequence = false;
ResIdSet selection;
char* selection_file = 0;
int size = 0;
bool verbose = false;
bool extract = false;
bool mcsym_file = false;

int main (int argc, char* argv[])
{
  gOut.setVerboseLevel (1);

  if ( argc == 1 )
    {
      cerr << "Usage: " << argv[0] << " ([-a] | [-e selection -f file_sel -s size] | [-m]) pdbfile+" << endl;
      exit (EXIT_SUCCESS);
    }

  while ((__current_option = getopt (argc, argv, __shortopts)) != EOF)
    {
      switch (__current_option)
	{
	case 'a':
	  __sequence = true;
	  break;
	case 'e':
	  selection = ResIdSet (optarg);
	  cerr << "selection = " << selection << endl;
	  extract = true;
	  break;
	case 'f':
	  selection_file = optarg;
	  extract = true;
	  break;
	case 'm':
	  mcsym_file = true;
	  break;
	case 's':
	  size = atoi (optarg);
	  break;	
	case 'v':
	  verbose = true;
	  break;  	  
	default:
	  cerr << "Usage: " << argv[0] << " [-e selection] pdbfile+" << endl;
	  exit (EXIT_FAILURE);
	}
    }
  
  
  clock_t startTime, endTime;
  startTime = clock();
  
  
  while (optind<argc) {
    
//     {
//       izfPdbstream fin (argv[optind]);
      
//       if (fin.good ()) {
// 	MolecularComplex *m = new MolecularComplex;
// 	fin >> *m;
	
// 	cout << *m << endl;
//       }
//       fin.close ();
//     }

    //xit (0);
    
    izfPdbstream fin (argv[optind]);
    PdbFileHeader header;
    Model model, model_orig;
    
    if (verbose) cout << argv[optind];

    if (!fin) {
      cerr << "File could not be opened" << endl;
      exit (EXIT_FAILURE); 
    }
    header = fin.getHeader ();
    fin >> model;
    model.removeWater ();
    fin.close ();
    
    if (verbose) cout << " (" << model.size () << " residues) " << flush;
    
    
    
    if (model.size () > 0) {    	      
      AnnotatedModel amodel (&model);      
      if (verbose) cerr << "Annotating..." << flush;
      amodel.annotate ();      
      
      
      if (mcsym_file != true) 
	{
	  amodel.dumpMcc (header.getPdbId ().c_str ());
	} 
      else 
	{	  	  
	  if (verbose) cerr << "Finding helices..." << flush;
	  amodel.findHelices ();
	  if (verbose) cerr << "done." << endl;
	  if (verbose) cerr << "Finding strands..." << flush;
	  amodel.findStrands ();
	  amodel.classifyStrands ();
	  
	  cout << "Residue conformations -------------------------------------------" << endl;      
	  amodel.dumpConformations ();
	  cout << endl;      
	  
	  amodel.dumpStacks ();
	  cout << endl;
	  
	  cout << "Base-pairs ------------------------------------------------------" << endl;	    
	  amodel.findKissingHairpins ();
	  cout << endl;     
	  
	  cout << "Triples ---------------------------------------------------------" << endl;
	  amodel.dumpTriples ();      
	  cout << endl;
	  
	  cout << "Helices ---------------------------------------------------------" << endl;
	  amodel.dumpHelices ();
	  cout << endl;
	  
	  cout << "Strands ---------------------------------------------------------" << endl;
	  amodel.dumpStrands ();
	  cout << endl;
	  
	  cout << "Various features ------------------------------------------------" << endl;
	  
	  if (verbose) cerr << "Finding pseudoknots..." << endl;
	  amodel.findPseudoknots ();
	  cout << endl;
	  
	  cout << "Sequences -------------------------------------------------------" << endl;
	  amodel.dumpSequences ();
	  cout << endl;	 	  
	}
    }
    
    optind++;
  }
  
  return EXIT_SUCCESS;
}


// 	if (extract) {
// 	  Model model_sel, model_sel_orig;
// 	  if (verbose) cout << "Extracting " << size << endl;

// 	  if (selection_file) {
// 	    izfPdbstream fin (selection_file);
// 	    if (!fin) {
// 	      cerr << "File could not be opened" << endl;
// 	    }
// 	    Model::iterator i;
// 	    fin >> model_sel;
// 	    fin.close ();
	    
// 	    model_sel_orig = model_sel;
// 	    model_sel.keepNucleicAcid ();
// 	    model_sel.validate ();
	    
// 	    for (i=model_sel.begin (); i!=model_sel.end (); ++i) {
// 	      selection.insert ((CResId&)*i);
// 	    }
// 	  }
// 	  CResIdSet ext_selection;
// 	  ext_selection = amodel.extract (selection, size);
	  
// 	  CResIdSet::iterator i;
// 	  Model::iterator j;
// 	  oPdbstream fout (cout.rdbuf ());
// 	  fout << "REMARK Selection: " << selection << endl;
// 	  fout << "REMARK Extended selection: " << ext_selection << endl;
// // 	  fout << "REMARK Original residues: " << endl;
// // 	  fout << model_sel_orig;
// 	  fout << "REMARK Extracted from " << argv[optind] << ":" << endl;
// 	  for (i=ext_selection.begin (); i!=ext_selection.end (); ++i) {
// 	    j = model_orig.find (*i);
// 	    if (j!=model_orig.end ()) {
// 	      fout << *j;
// 	    }
// 	  }
// 	  fout << "END\n";	 
// 	} else {






// 	  if (__sequence) {
// 	    cout << argv[optind] << " : " << flush;
// 	    amodel.dumpSequences (false);
// 	    cout << endl;
// 	    amodel.dumpPairs ();
// 	  }
// 	  else {
// //  	    cerr << "Block decomposition..." << endl;
	    
// //   	    amodel.decompose ();
// // //  	    amodel.dumpSequences ();
	    
// //   	    amodel.dumpPDF ("out.pdf");

// 	    //	    amodel.dumpSimplePDF (argv[optind], "out.pdf");
	    	     
// 	  }
// 	}
//       }
//       endTime = clock();
//       //     cerr << "Elapsed time: " << ( (float)(endTime - startTime) / (float)CLOCKS_PER_SEC ) << " sec." << endl;
      
//     }

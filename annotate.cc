//                              -*- Mode: C++ -*- 
// annote.cc
// Copyright © 2001, 2002, 2003, 2004 Laboratoire de Biologie Informatique et Théorique.
// Author           : Patrick Gendron
// Created On       : Fri May 18 09:38:07 2001
// Last Modified By : Philippe Thibault
// Last Modified On : Tue Jul  6 14:57:34 2004
// Update Count     : 128
// Status           : Unknown.
// 

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream>
#include <stdlib.h>
#include <malloc.h>
#include <getopt.h>

#include "annotate.h"

#include "mccore/Pdbstream.h"
#include "mccore/Binstream.h"
#include "mccore/Residue.h"
#include "mccore/Model.h"
#include "mccore/Algo.h"
#include "mccore/Vector3D.h"
#include "mccore/Messagestream.h"
#include "mccore/zstream.h"
#include "mccore/Relation.h"
#include "mccore/PropertyType.h"
#include "mccore/stlio.h"
#include "mccore/Exception.h"

#include "AnnotatedModel.h"

#include "MolecularComplex.h"

#define NONE_MARK  0
#define HELIX_MARK 1
#define PAIR_MARK  2
#define BULGE_MARK 3
#define LOOP_MARK  4

int __current_option;
const char* __shortopts = "e:f:hHms:vV?";
ResIdSet selection;
char* selection_file = 0;
int size = 0;
bool verbose = false;
bool extract = false;
bool mcsym_file = false;

/**
 * @short Version function
 * @author Patrick Gendron
 * -----------------------------------------------------------------------------
 */
void __version (void)
{
  cout << prog_name << " " << prog_major_version << "." << prog_minor_version
       << " (" << __DATE__ << ")" << endl
       << "  using lib" << "mccore" << " " << VERSION << endl;
}


/**
 * @short Usage function
 * @author Patrick Gendron
 * -----------------------------------------------------------------------------
 */
void __usage (void)
{
  cout << "usage: " << prog_name 
       << " [-hHmvV?] [-e selection -f file_sel -s size] [structure file]" << endl;
}

/**
 * @short Help function
 * @author Patrick Gendron
 * -----------------------------------------------------------------------------
 */
void __help (void)
{
  cout << "This program annotates a structure (and more)." << endl
       << "  -e sel  extract these residues from the structure" << endl 
       << "  -f pdb  extract residues found in this pdb file from the structure" << endl 
       << "  -h      help      print this help" << endl
       << "  -m      produce an output in McSym format (with backtrack)" << endl
       << "  -s N    extend the selection stated with -e or -f by N connected residues in the graph" << endl
       << "  -v      verbose   be verbose" << endl
       << "  -V      version   print the software version info" << endl;
}




/**
 * @short Annotate main algorithm.
 * @author Patrick Gendron
 * -----------------------------------------------------------------------------
 */
int main (int argc, char* argv[])
{
  gOut.setVerboseLevel (1);

  if ( argc == 1 )
  {
    __usage ();
    exit (EXIT_SUCCESS);
  }

  while ((__current_option = getopt (argc, argv, __shortopts)) != EOF)
  {
    switch (__current_option)
    {
    case 'e':
      selection = ResIdSet (optarg);
      cerr << "selection = " << selection << endl;
      extract = true;
      break;
    case 'f':
      selection_file = optarg;
      extract = true;
      break;
    case 'h': case 'H':
      __usage (); __help (); exit (EXIT_SUCCESS); break;
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
  
  try
  {
  
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
      model.addHLP ();
      fin.close ();
    
      if (verbose) cout << " (" << model.size () << " residues) " << flush;
    
    
    
      if (model.size () > 0) {    	      
	AnnotatedModel amodel (&model);      
	if (verbose) cerr << "Annotating..." << flush;
	amodel.annotate ();      
      
      
	if (mcsym_file == true) 
	{
	  amodel.dumpMcc (header.getPdbId ().c_str ());
	} 
	else if (extract) 
	{
	  Model model_sel, model_sel_orig;
	  if (verbose) cout << "Extracting " << size << endl;
	
	  if (selection_file) {
	    izfPdbstream fin (selection_file);
	    if (!fin) {
	      cerr << "File could not be opened" << endl;
	    }
	    Model::iterator i;
	    fin >> model_sel;
	    fin.close ();
	  
	    model_sel_orig = model_sel;
	    model_sel.keepNucleicAcid ();
	    model_sel.validate ();
	  
	    for (i=model_sel.begin (); i!=model_sel.end (); ++i) {
	      selection.insert (i->getResId ());
	    }
	  }
	  ResIdSet ext_selection;
	  ext_selection = amodel.extract (selection, size);
	
	  ResIdSet::iterator i;
	  Model::iterator j;
	  oPdbstream fout (cout.rdbuf ());
	  fout << "REMARK Selection: " << selection << endl;
	  fout << "REMARK Extended selection: " << ext_selection << endl;
	  // 	  fout << "REMARK Original residues: " << endl;
	  // 	  fout << model_sel_orig;
	  fout << "REMARK Extracted from " << argv[optind] << ":" << endl;
	  for (i=ext_selection.begin (); i!=ext_selection.end (); ++i) {
	    j = model.find (*i);
	    if (j!=model.end ()) {
	      fout << *j;
	    }
	  }
	  fout << "END\n";	 
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

  }
  catch (Exception ex)
  {
    cerr << argv[0] << ": " << ex << endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

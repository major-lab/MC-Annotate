//                              -*- Mode: C++ -*- 
// annotate.cc
// Copyright © 2001-04 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri May 18 09:38:07 2001
// $Revision$
// $Id$


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream>
#include <cstdlib>
#include <string>
#include <unistd.h>

//#include "mccore/Algo.h"
//#include "mccore/Binstream.h"
#include "mccore/Exception.h"
#include "mccore/Messagestream.h"
#include "mccore/GraphModel.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
//#include "mccore/Relation.h"
//#include "mccore/Residue.h"
#include "mccore/RnamlReader.h"
//#include "mccore/PropertyType.h"
//#include "mccore/Vector3D.h"
//#include "mccore/stlio.h"
//#include "mccore/zstream.h"

#include "AnnotateModule.h"
//#include "AnnotatedModel.h"
//#include "MolecularComplex.h"

using namespace mcannotate;
using namespace mccore;
using namespace std;

#define NONE_MARK  0
#define HELIX_MARK 1
#define PAIR_MARK  2
#define BULGE_MARK 3
#define LOOP_MARK  4

const char *shortopts = "Ve:f:hms:v";
ResIdSet selection;
string selection_file;
int size = 0;
bool mcsym_file = false;



void
version ()
{
  gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << endl;
}


void
usage ()
{
  gOut (0) << "usage: " << PACKAGE
	   << " [-hmvV] [-e selection -f file_sel -s size] [structure file]" << endl;
}


void
help ()
{
  gOut (0) << "This program annotates a structure (and more)." << endl
	   << "  -e sel  extract these residues from the structure" << endl 
	   << "  -f pdb  extract residues found in this pdb file from the structure" << endl 
	   << "  -h      help      print this help" << endl
	   << "  -m      produce an output in McSym format (with backtrack)" << endl
	   << "  -s N    extend the selection stated with -e or -f by N connected residues in the graph" << endl
	   << "  -v      verbose   be verbose" << endl
	   << "  -V      version   print the software version info" << endl;
}


void
read_options (int argc, char *argv[])
{
  int c;
  
  while ((c = getopt (argc, argv, shortopts)) != EOF)
    {
      switch (c)
	{
	case 'V':
	  version ();
	  exit (EXIT_SUCCESS);
	  break;
	case 'e':
	  selection = ResIdSet (optarg);
	  break;
	case 'f':
	  selection_file = optarg;
	  break;
	case 'h':
	  usage ();
	  help ();
	  exit (EXIT_SUCCESS);
	  break;
	case 'm':
	  mcsym_file = true;
	  break;
	case 's':
	  size = atoi (optarg);
	  break;	
	case 'v':
	  gOut.setVerboseLevel (gOut.getVerboseLevel () + 1);
	  break;
	default:
	  usage ();
	  exit (EXIT_FAILURE);
	}
    }
  if (argc <= optarg)
    {
      usage ();
      exit (EXIT_FAILURE);
    }
}


void
addFile (AnnotateModule *module, const char *name)
{
  RnamlReader reader (name);
  Molecule *molecule;
  
  if (0 == (molecule = reader.read ()))
    {
      izfPdbstream fin (name);
      PdbFileHeader header;

      if (! fin)
	{
	  cerr << "File could not be opened" << endl;
	  exit (EXIT_FAILURE); 
	}
      header = fin.getHeader ();
      while (!fin)
	{
	  GraphModel model;
	  
	  fin >> model;
	  module->add (name, header, model);
	}
      fin.close ();
    }
  else
    {
      Molecule::iterator it;
      
      for (it = molecule->begin (); molecule->end () != it; ++it, ++count)
	{
	  module->add (name, header, *it);
	}
      delete molecule;
    }
  reader.close ();
}
  

int
main (int argc, char* argv[])
{
  read_options (argc, argv);

  try
    {
      AnnotateModule *module;
      View *view;
      
      module = new AnnotateModule (selection, selection_file, size);
      
      while (optind < argc)
	{
	  addFile (module, argv[optind++]);
	}
      module->exec ();

      if (mcsym_file)
	{
	  module->writeMCSym (cout);
	}
      if (! rnamlOutput.empty ())
	{
	  module->writeRNAML (rnamlOutput);
	}
      view = new TerminalView (module, gOut);
      view->show ();

      delete view;
      delete module;
    }
  catch (Exception& ex)
    {
      gOut (1) << PACKAGE << ": " << ex << endl;
      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
}

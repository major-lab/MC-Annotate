//                              -*- Mode: C++ -*- 
// annotate.cc
// Copyright © 2001-05 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri May 18 09:38:07 2001
// $Revision$
// $Id$


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <iostream>
#include <string>
#include <unistd.h>

#include "mccore/Messagestream.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"

#include "AnnotateModel.h"

#undef PACKAGE
#undef VERSION
#define PACKAGE "annotate"
#define VERSION "1.0-mccore1.5"

using namespace mccore;
using namespace std;
using namespace annotate;

#define NONE_MARK  0
#define HELIX_MARK 1
#define PAIR_MARK  2
#define BULGE_MARK 3
#define LOOP_MARK  4

const char* shortopts = "e:f:hlms:vV";
ResIdSet selection;
string selection_file;
unsigned int size = 0;
bool verbose = false;
bool extract = false;
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
    << " [-hlmvV] [-e selection -f file_sel -s size] [structure file]"
    << endl;
}


void
help ()
{
  gOut (0)
    << "This program annotates a structure (and more)." << endl
    << "  -e sel            extract these residues from the structure" << endl 
    << "  -f pdb            extract residues found in this pdb file from the structure" << endl 
    << "  -h      help      print this help" << endl
    << "  -l      verbose   be more verbose (log)" << endl
    << "  -m                produce an output in McSym format (with backtrack)" << endl
    << "  -s N              extend the selection stated with -e or -f by N connected residues in the graph" << endl
    << "  -v      verbose   be verbose" << endl
    << "  -V      version   print the software version info" << endl;
}


void
read_options (int argc, char* argv[])
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
          gErr (3) << "selection = " << selection << endl;
          extract = true;
          break;
        case 'f':
          selection_file = optarg;
          extract = true;
          break;
        case 'h':
          usage ();
          help ();
          exit (EXIT_SUCCESS);
          break;
        case 'l':
          gErr.setVerboseLevel (gErr.getVerboseLevel () + 1);
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

  if (argc - optind < 1)
    {
      usage ();
      exit (EXIT_FAILURE);
    }
}



int
main (int argc, char *argv[])
{
  read_options (argc, argv);

  while (optind < argc)
    {
      izfPdbstream is;

      is.open (argv[optind]);
      if (! is.good ())
	{
	  gErr (0) << PACKAGE << ": cannot open pdb file '"
		   << argv[optind] << "'." << endl;
	}
      else
	{
	  AnnotateModelFM afm;
	  Molecule molecule (&afm);
	  Molecule::iterator molIt;

	  gErr (3) << PACKAGE << ": reading " << argv[optind] << endl;   
	  is >> molecule;
	  is.close ();
	  for (molIt = molecule.begin (); molecule.end () != molIt; ++molIt)
	    {
	      AnnotateModel &am = (AnnotateModel&) *molIt;

	      am.annotate ();
	      if (mcsym_file == true) 
		{
		  // NEEDS TO BE FIXED!
		  //            PdbFileHeader header;
		  //            header = is.getHeader ();
		  am.dumpMcc (" "); //header.getPdbId ().c_str ()
		} 
	      else if (extract) 
		{
		  // NEEDS TO BE FIXED!
		}
	      else
		{
		  gOut(0) << am ;
		}
	    }
	}
      ++optind;
    }
  return EXIT_SUCCESS;	
}

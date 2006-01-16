//                              -*- Mode: C++ -*- 
// annotate.cc
// Copyright © 2001-06 Laboratoire de Biologie Informatique et Théorique.
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

#include "mccore/Binstream.h"
#include "mccore/Messagestream.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResidueFactoryMethod.h"
#ifdef HAVE_LIBRNAMLC__
#include "mccore/RnamlReader.h"
#endif

#include "AnnotateModel.h"

using namespace mccore;
using namespace std;
using namespace annotate;

bool binary = false;
bool extract = false;
bool mcsym_file = false;
ResIdSet selection;
string selection_file;
const char* shortopts = "Vbe:f:hlms:v";
unsigned int size = 0;



void
version ()
{
  gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << endl;
}


void
usage ()
{
  gOut (0) << "usage: " << PACKAGE
    << " [-bhlmvV] [-e selection -f file_sel -s size] [structure file]"
    << endl;
}


void
help ()
{
  gOut (0)
    << "This program annotates a structure (and more)." << endl
    << "  -b                read binary files instead of pdb files" << endl
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
	case 'b':
	  binary = true;
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


mccore::Molecule*
loadFile (const string &filename)
{
  Molecule *molecule;

  molecule = 0;
  if (binary)
    {
      izfBinstream in;

      in.open (filename.c_str ());
      if (in.fail ())
	{
	  gErr (0) << PACKAGE << ": cannot open binary file '" << filename << "'." << endl;
	  return 0;
	}
      molecule = new Molecule ();
      in >> *molecule;
      in.close ();
    }
  else
    {
#ifdef HAVE_LIBRNAMLC__
      RnamlReader reader (filename.c_str ());
      
      if (0 == (molecule = reader.read ()))
#endif
	{
	  izfPdbstream in;
	  ResidueFM rFM;
	  GraphModelFM gFM (&rFM);
	  
	  in.open (filename.c_str ());
	  if (in.fail ())
	    {
	      gErr (0) << PACKAGE << ": cannot open pdb file '" << filename << "'." << endl;
	      return 0;
	    }
	  molecule = new Molecule (&gFM);
	  in >> *molecule;
	  in.close ();
	}
    }
  return molecule;
}


int
main (int argc, char *argv[])
{
  read_options (argc, argv);

  while (optind < argc)
    {
      Molecule *molecule;
      Molecule::iterator molIt;
      
      gErr (3) << PACKAGE << ": reading " << argv[optind] << endl;
      molecule = loadFile (argv[optind]);
      if (0 != molecule)
	{
	  for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
	    {
	      AnnotateModel am (*molIt);

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
	  delete molecule;
	}
      ++optind;
    }
  return EXIT_SUCCESS;	
}

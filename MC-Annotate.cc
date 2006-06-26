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
#include "mccore/Version.h"

#include "AnnotateModel.h"

using namespace mccore;
using namespace std;
using namespace annotate;

bool binary = false;
unsigned int environment = 0;
bool oneModel = false;
unsigned int modelNumber = 0;  // 1 based vector identifyer, 0 means all
string residueSelection = "";
#ifdef HAVE_LIBRNAMLC__
const char* shortopts = "Vbe:f:hlr:vx";
bool xmlFile = false;
#else
const char* shortopts = "Vbe:f:hlr:v";
#endif



void
version ()
{
  mccore::Version mccorev;
  
  gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << endl
	   << "  using " << mccorev << endl;
}


void
usage ()
{
  gOut (0) << "usage: " << PACKAGE
	   << " [-bhlvV] [-e num] [-f <model number>] [-r <residue ids>]"
#ifdef HAVE_LIBRNAMLC__
	   << " [-x]"
#endif
	   << " structure file ..." << endl;
}


void
help ()
{
  gOut (0)
    << "This program annotates a structure (and more)." << endl
    << "  -b                read binary files instead of pdb files" << endl
    << "  -e num            number of surrounding layers of connected residues to annotate" << endl
    << "  -f model number   model to print" << endl
    << "  -h                print this help" << endl
    << "  -l                be more verbose (log)" << endl
    << "  -r sel            extract these residues from the structure" << endl 
    << "  -v                be verbose" << endl
    << "  -V                print the software version info" << endl
#ifdef HAVE_LIBRNAMLC__
    << "  -x                save the annotation in rnaml format" << endl
#endif
    ;    
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
	  {
	    int tmp;

	    tmp = atoi (optarg);
	    if (0 < tmp)
	      {
		environment = tmp;
	      }
	  }
	  break;
	case 'f':
	  {
	    int tmp;
	   
	    tmp = atoi (optarg);
	    if (0 <= tmp)
	      {
		modelNumber = tmp;
		oneModel = true;
	      }
	  }
	  break;
        case 'h':
          usage ();
          help ();
          exit (EXIT_SUCCESS);
          break;
        case 'l':
          gErr.setVerboseLevel (gErr.getVerboseLevel () + 1);
          break;
	case 'r':
	  residueSelection = optarg;
	  break;
	case 'v':
	  gOut.setVerboseLevel (gOut.getVerboseLevel () + 1);
          break;
#ifdef HAVE_LIBRNAMLC__
	case 'x':
	  xmlFile = true;
	  break;
#endif
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
  ResidueFM rFM;
  AnnotateModelFM aFM (residueSelection, environment, &rFM);

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
      molecule = new Molecule (&aFM);
      in >> *molecule;
      in.close ();
    }
  else
    {
#ifdef HAVE_LIBRNAMLC__
      RnamlReader reader (filename.c_str (), &aFM);
      
      if (0 == (molecule = reader.read ()))
	{
#endif
	  izfPdbstream in;
	  
	  in.open (filename.c_str ());
	  if (in.fail ())
	    {
	      gErr (0) << PACKAGE << ": cannot open pdb file '" << filename << "'." << endl;
	      return 0;
	    }
	  molecule = new Molecule (&aFM);
	  in >> *molecule;
	  in.close ();
#ifdef HAVE_LIBRNAMLC__
	}
#endif
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
      
      molecule = loadFile ((string) argv[optind]);
      if (0 != molecule)
	{
	  for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
	    {
	      if (0 != modelNumber)
		{
		  --modelNumber;
		}
	      else
		{
		  AnnotateModel &am = (AnnotateModel&) *molIt;
		  
		  am.annotate (true);
		  gOut(0) << am ;
		  if (oneModel)
		    {
		      optind = argc - 1;
		      break;
		    }
		}
	    }
	  delete molecule;
	}
      ++optind;
    }
  return EXIT_SUCCESS;	
}

//                              -*- Mode: C++ -*- 
// annotate.cc
// Copyright © 2001-06 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Patrick Gendron
// Created On       : Fri May 18 09:38:07 2001
// $Revision: 58 $
// $Id: MC-Annotate.cc 58 2006-11-15 21:09:19Z larosem $


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cerrno>
#include <cstdlib>
#include <string>
#include <sstream>
#include <unistd.h>

#include "mccore/Binstream.h"
#include "mccore/Exception.h"
#include "mccore/Messagestream.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/Molecule.h"
#include "mccore/Pdbstream.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResidueFactoryMethod.h"
#include "mccore/ResIdSet.h"
#ifdef HAVE_LIBRNAMLC__
#include "mccore/RnamlReader.h"
#endif
#include "mccore/Version.h"

#include "AnnotateModel.h"
#include "AnnotationCycles.h"
#include "AnnotationLinkers.h"
#include "AnnotationLoops.h"
#include "AnnotationStems.h"
#include "AnnotationTertiaryCycles.h"
#include "AnnotationTertiaryPairs.h"

using namespace mccore;
using namespace std;
using namespace annotate;

bool binary = false;
unsigned int environment = 0;
bool oneModel = false;
unsigned int modelNumber = 0;  // 1 based vector identifier, 0 means all
std::string gstrOutputDirectory = "";
ResIdSet residueSelection;
const char* shortopts = "Vbd:e:f:hlr:v";



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
	   << " [-bhlvV] [-e num] [-d <output directory>] [-f <model number>] [-r <residue ids>] <structure file> ..."
	   << endl;
}


void
help ()
{
  gOut (0)
    << "This program annotate structures (and more)." << endl
    << "  -b                read binary files instead of pdb files" << endl
    << "  -d <directory>    output directory" << endl
    << "  -e num            number of surrounding layers of connected residues to annotate" << endl
    << "  -f model number   model to print" << endl
    << "  -h                print this help" << endl
    << "  -l                be more verbose (log)" << endl
    << "  -r sel            extract these residues from the structure" << endl 
    << "  -v                be verbose" << endl
    << "  -V                print the software version info" << endl;    
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
	    long int tmp;

	    tmp = strtol (optarg, 0, 10);
	    if (ERANGE == errno
		|| EINVAL == errno
		|| 0 > tmp)
	      {
		gErr (0) << PACKAGE << ": invalid environment value." << endl;
		exit (EXIT_FAILURE);
	      }
	    environment = tmp;
	    break;
	  }
	case 'd':
	{
		gstrOutputDirectory = optarg;
		if (ERANGE == errno	|| EINVAL == errno)
		{
			gErr (0) << PACKAGE << ": invalid environment value." << endl;
			exit (EXIT_FAILURE);
		}
		break;		
	}
	case 'f':
	  {
	    long int tmp;

	    tmp = strtol (optarg, 0, 10);
	    if (ERANGE == errno
		|| EINVAL == errno
		|| 0 > tmp)
	      {
		gErr (0) << PACKAGE << ": invalid model value." << endl;
		exit (EXIT_FAILURE);
	      }
	    modelNumber = tmp;
	    oneModel = true;
	    break;
	  }
        case 'h':
          usage ();
          help ();
          exit (EXIT_SUCCESS);
          break;
        case 'l':
          gErr.setVerboseLevel (gErr.getVerboseLevel () + 1);
          break;
	case 'r':
	  try
	    {
	      residueSelection.insert (optarg);
	    }
	  catch (IntLibException &e)
	    {
	      gErr (0) << PACKAGE << ": invalid residue selection." << endl;
	      exit (EXIT_FAILURE);
	    }
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

std::string getFilePrefix(const std::string& aFileName)
{
	std::string::size_type index;
	std::string filename = aFileName;
	if (std::string::npos != (index = filename.rfind ("/")))
    {
		filename.erase (0, index + 1);
    }
	if (string::npos != (index = filename.find (".")))
    { 
		filename.erase (index, filename.size ());
    }
	return filename;
}

void dumpCyclesFiles(
	const std::string& aFilePrefix, 
	const AnnotationTertiaryCycles& aAnnotationCycles)
{
	int i = 1;
	for(std::vector<annotate::Cycle>::const_iterator it = aAnnotationCycles.getCycles().begin();
		it != aAnnotationCycles.getCycles().end();
		++ it)
	{
		std::ostringstream oss;
		if(!gstrOutputDirectory.empty())
		{
			oss << gstrOutputDirectory << "/";
		}
		oss << aFilePrefix << "_c" << i << ".pdb.gz";
		ozfPdbstream ops;
		ops.open ((oss.str ()).c_str ());
		if (! ops)
		{
			gErr (0) << "Cannot open file " << oss.str() << endl;
			continue;
		}
		ops << *it;
		ops.close ();
		++ i;
	}	
}


int
main (int argc, char *argv[])
{
	read_options (argc, argv);

	while (optind < argc)
    {
		Molecule *molecule;
		Molecule::iterator molIt;
		std::string filename = (std::string) argv[optind];
      
		molecule = loadFile (filename);
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
					AnnotationStems annStems;
					AnnotationLinkers annLinkers;
					AnnotationLoops annLoops;
					AnnotationTertiaryPairs annTertiaryPairs;
					AnnotationCycles annCycles;
					AnnotationTertiaryCycles annTertiaryCycles;
		  
					am.addAnnotation(annStems);
					am.addAnnotation(annLinkers);
					am.addAnnotation(annLoops);
					am.addAnnotation(annTertiaryPairs);
					am.addAnnotation(annCycles);
					am.addAnnotation(annTertiaryCycles);
					am.annotate ();
					
					gOut (0) << "Annotating Model ------------------------------------------------" << endl;
					gOut (0) << filename << std::endl;
					gOut(0) << am;
					
					dumpCyclesFiles(getFilePrefix(filename), annTertiaryCycles);					
					
					if (oneModel)
					{
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

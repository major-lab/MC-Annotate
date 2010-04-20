//                              -*- Mode: C++ -*-
// helix_angle.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Mon Mar 22 10:22:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "helix_angle.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>

static const char* shortopts = "Vhlvc:p:";

HelixAngle::HelixAngle(int argc, char * argv [])
{
	// Read the command line options
	readOptions(argc, argv);


}

void HelixAngle::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void HelixAngle::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-hlvV] <structure file> ..."
		<< std::endl;
}

void HelixAngle::help () const
{
	mccore::gOut (0)
		<< "This program computes the phi-phi angles between cyclinders associated with a fixed and a mobile stem" << std::endl
		<< "(Reproducing the mesure presented in Chu and al. 2009 RNA." << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
}


void HelixAngle::readOptions (int argc, char* argv[])
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
		case 'h':
			usage ();
			help ();
			exit (EXIT_SUCCESS);
			break;
		case 'l':
			mccore::gErr.setVerboseLevel (mccore::gErr.getVerboseLevel () + 1);
			break;
		case 'v':
			mccore::gOut.setVerboseLevel (mccore::gOut.getVerboseLevel () + 1);
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

mccore::Molecule* HelixAngle::loadFile (const std::string &astrFilename)
{
  Molecule *molecule;
  ResidueFM rFM;
  annotate::AnnotateModelFM aFM (residueSelection, environment, &rFM);

  molecule = 0;
  if (binary)
    {
      izfBinstream in;

      in.open (astrFilename.c_str ());
      if (in.fail ())
	{
	  gErr (0) << PACKAGE << ": cannot open binary file '" << astrFilename << "'." << endl;
	  return 0;
	}
      molecule = new Molecule (&aFM);
      in >> *molecule;
      in.close ();
    }
  else
    {
#ifdef HAVE_LIBRNAMLC__
      RnamlReader reader (astrFilename.c_str (), &aFM);

      if (0 == (molecule = reader.read ()))
	{
#endif
	  izfPdbstream in;

	  in.open (astrFilename.c_str ());
	  if (in.fail ())
	    {
	      gErr (0) << PACKAGE << ": cannot open pdb file '" << astrFilename << "'." << endl;
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

void HelixAngle::processModels() const
{
	while (optind < argc)
	{
		Molecule *molecule;
		Molecule::iterator molIt;
		std::string filename = (std::string) argv[optind];

		molecule = loadFile (filename);
		if (0 != molecule)
		{
			unsigned int uiCurrentModel = 1;
			for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
			{
				if (0 != modelNumber)
				{
					--modelNumber;
				}
				else
				{
					annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
					std::string strPDBPrefix = getFilePrefix(filename);
					am.name(strPDBPrefix);
					am.id(uiCurrentModel);

					if()
					{
						std::cout << strPDBPrefix << " : ";
						std::cout << uiCurrentModel << " : ";
						std::cout << "(" << fFixFace;
						std::cout << "," << fMobileFace << ")";
						std::cout << std::endl;
					}

//					am.annotate (gucRelationMask);
					if (oneModel)
					{
						break;
					}
				}
				uiCurrentModel ++;
			}
			delete molecule;
		}
		++optind;
	}
}

int main (int argc, char *argv[])
{
	HelixAngle theApp(argc, argv);

	return EXIT_SUCCESS;
}

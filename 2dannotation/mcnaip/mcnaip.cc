//                              -*- Mode: C++ -*-
// mcnaip.cc
// Copyright © 2009 Laboratoire de Biologie Informatique et Théorique.
//                  Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Aug 18 11:38:07 2009
// $Revision: 1 $
// $Id: mcnaip.cc 1 2009-09-18 11:38:07 blanchmf $



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
#include "AnnotationChains.h"
#include "AnnotationInteractions.h"
#include "AnnotationLinkers.h"
#include "AnnotationLoops.h"
#include "AnnotationResSecondaryStructures.h"
#include "AnnotationStems.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotationTertiaryStacks.h"
#include "StringTable.h"

bool binary = false;
unsigned int environment = 0;
bool oneModel = false;
unsigned int guiModelNumber = 0;  // 1 based vector identifier, 0 means all
unsigned char gucRelationMask =
	mccore::Relation::adjacent_mask
	| mccore::Relation::pairing_mask
	| mccore::Relation::stacking_mask
	| mccore::Relation::backbone_mask;
const char* shortopts = "Vbf:hlvm:";

void version ()
{
	mccore::Version mccorev;

	mccore::gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
	   << "  using " << mccorev << std::endl;
}


void usage ()
{
	mccore::gOut (0) << "usage: " << PACKAGE
	   << " [-bhlvV] [-e num] [-f <model number>] [-m <relation mask>] <structure file> ..."
	   << std::endl;
}

void help ()
{
	mccore::gOut(0)
		<< "This program annotate structures, compute a 2D layount and returns " << std::endl
		<< "the list of interactions between non adjacent secondary structure " << std::endl
		<< "components." << std::endl
		<< "  -b                read binary files instead of pdb files" << std::endl
		<< "  -m <mask>         annotation mask string: any combination of " << std::endl
		<< "'A' (adjacent), 'S' (stacking), 'P' (pairing) and 'B' (backbone). (default: all)" << std::endl
		<< "  -f model number   model to print" << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
}


void read_options (int argc, char* argv[])
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
			if (ERANGE == errno || EINVAL == errno || 0 > tmp)
			{
				mccore::gErr (0) << PACKAGE << ": invalid environment value." << std::endl;
				exit (EXIT_FAILURE);
			}
			environment = tmp;
			break;
		}
		case 'f':
		{
			long int tmp;

			tmp = strtol (optarg, 0, 10);
			if (ERANGE == errno || EINVAL == errno || 0 > tmp)
			{
				mccore::gErr (0) << PACKAGE << ": invalid model value." << std::endl;
				exit (EXIT_FAILURE);
			}
			guiModelNumber = tmp;
			oneModel = true;
			break;
		}
		case 'h':
			usage ();
			help ();
			exit (EXIT_SUCCESS);
			break;
		case 'l':
			mccore::gErr.setVerboseLevel (mccore::gErr.getVerboseLevel () + 1);
			break;
		case 'm':
		{
			std::string strMask = optarg;
			gucRelationMask = 0;
			if(strMask.find("A") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::adjacent_mask;
    		}
    		if(strMask.find("S") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::stacking_mask;
    		}
    		if(strMask.find("P") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::pairing_mask;
    		}
    		if(strMask.find("B") != std::string::npos)
    		{
    			gucRelationMask |= mccore::Relation::backbone_mask;
    		}
			break;
    	}
		case 'v':
			mccore::gOut.setVerboseLevel (gOut.getVerboseLevel () + 1);
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


mccore::Molecule* loadFile (const std::string &filename)
{
	mccore::ResIdSet residueSelection;
	mccore::Molecule *molecule;
	mccore::ResidueFM rFM;
	annotate::AnnotateModelFM aFM (residueSelection, environment, &rFM);

	molecule = 0;
	if (binary)
	{
		mccore::izfBinstream in;

		in.open (filename.c_str());
		if (in.fail ())
		{
			mccore::gErr(0) << PACKAGE << ": cannot open binary file '";
			mccore::gErr(0) << filename << "'." << endl;
			return 0;
		}
		molecule = new mccore::Molecule (&aFM);
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
				mccore::gErr(0) << PACKAGE << ": cannot open pdb file '";
				mccore::gErr(0) << filename << "'." << std::endl;
				return 0;
			}
			molecule = new mccore::Molecule (&aFM);
			in >> *molecule;
			in.close ();
#ifdef HAVE_LIBRNAMLC__
		}
#endif
	}
	return molecule;
}

std::string getPdbFileName(const std::string& aFileName)
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

std::string getModelString(unsigned int auiModel)
{
	std::ostringstream oss;
	oss << auiModel;
	return oss.str();
}

std::string getPairString(const annotate::BasePair& aPair)
{
	std::ostringstream oss;
	oss << aPair.fResId << "-" << aPair.rResId;
	return oss.str();
}

int main (int argc, char *argv[])
{
	annotate::StringTable stringTable(3);
	read_options (argc, argv);

	while (optind < argc)
    {
		Molecule *molecule;
		Molecule::iterator molIt;
		std::string filename = (std::string) argv[optind];

		molecule = loadFile (filename);
		if (0 != molecule)
		{
			unsigned int uiModel = 1;
			for (molIt = molecule->begin (); molecule->end () != molIt; ++molIt)
			{
				if (guiModelNumber != 0 && uiModel != guiModelNumber)
				{
					++ uiModel;
				}
				else
				{
					annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
					am.name(getPdbFileName(filename));
					annotate::AnnotationInteractions annInteractions;
					annotate::AnnotationChains annChains;
					annotate::AnnotationStems annStems;
					annotate::AnnotationLinkers annLinkers;
					annotate::AnnotationLoops annLoops;
					annotate::AnnotationTertiaryPairs annTertiaryPairs;
					annotate::AnnotationTertiaryStacks annTertiaryStacks;

		  			am.addAnnotation(annInteractions);
		  			am.addAnnotation(annChains);
					am.addAnnotation(annStems);
					am.addAnnotation(annLinkers);
					am.addAnnotation(annLoops);
					am.addAnnotation(annTertiaryPairs);
					am.addAnnotation(annTertiaryStacks);

					am.annotate (gucRelationMask);
					std::set<annotate::BasePair>::const_iterator itPair;
					for(itPair = annTertiaryPairs.getPairs().begin();
						itPair != annTertiaryPairs.getPairs().end();
						++ itPair)
					{
						std::vector<std::string>& tableRow = stringTable.addRow();
						tableRow[0] = getPdbFileName(filename);
						tableRow[1] = getModelString(uiModel);
						tableRow[2] = getPairString(*itPair);
					}

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

	mccore::gOut(0) << stringTable.toString(":");

	return EXIT_SUCCESS;
}

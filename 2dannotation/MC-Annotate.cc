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

#include <mccore/Binstream.h>
#include <mccore/Exception.h>
#include <mccore/Messagestream.h>
#include <mccore/ModelFactoryMethod.h>
#include <mccore/Molecule.h>
#include <mccore/Pdbstream.h>
#include <mccore/PropertyType.h>
#include <mccore/Relation.h>
#include <mccore/ResidueFactoryMethod.h>
#include <mccore/ResIdSet.h>
#ifdef HAVE_LIBRNAMLC__
#include <mccore/RnamlReader.h>
#endif
#include <mccore/Version.h>

#include "AnnotateModel.h"
#include "AnnotatedModelSQLExporter.h"
#include "AnnotationChains.h"
#include "AnnotationInteractions.h"
#include "AnnotationLinkers.h"
#include "AnnotationLoops.h"
#include "AnnotationResSecondaryStructures.h"
#include "AnnotationStems.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotationTertiaryStacks.h"

bool binary = false;
unsigned int environment = 0;
bool oneModel = false;
unsigned int modelNumber = 0;  // 1 based vector identifier, 0 means all
unsigned char gucRelationMask =
	mccore::Relation::adjacent_mask
	| mccore::Relation::pairing_mask
	| mccore::Relation::stacking_mask
	| mccore::Relation::backbone_mask;
ResIdSet residueSelection;
std::string gstrSQLFile = "";
const char* shortopts = "Vbe:f:hlr:vm:s:";

void
version ()
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void
usage ()
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-bhlvV] [-e num] [-f <model number>] [-r <residue ids>] [-m <relation mask>] [-s <SQL file>] <structure file> ..."
		<< std::endl;
}

void
help ()
{
	mccore::gOut (0)
		<< "This program annotate structures (and more)." << std::endl
		<< "  -b                read binary files instead of pdb files" << std::endl
		<< "  -m <mask>         annotation mask string: any combination of 'A' (adjacent), 'S' (stacking), 'P' (pairing) and 'B' (backbone). (default: all)" << std::endl
		<< "  -e num            number of surrounding layers of connected residues to annotate" << std::endl
		<< "  -f model number   model to print" << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
		<< "  -r sel            extract these residues from the structure" << std::endl
		<< "  -s <SQL file> 	export the annotation in SQL as relational database" << std::endl
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
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
			if (ERANGE == errno || EINVAL == errno || 0 > tmp)
			{
				mccore::gErr(0) << PACKAGE << ": invalid environment value.";
				mccore::gErr(0) << std::endl;
				exit (EXIT_FAILURE);
			}
			environment = tmp;
			break;
		}
		case 'f':
		{
			long int tmp;

			tmp = strtol (optarg, 0, 10);
			if (ERANGE == errno	|| EINVAL == errno || 0 > tmp)
			{
				mccore::gErr(0) << PACKAGE << ": invalid model value.";
				mccore::gErr(0) << std::endl;
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
		case 'r':
			try
			{
				residueSelection.insert (optarg);
			}
			catch (IntLibException &e)
			{
				mccore::gErr(0) << PACKAGE << ": invalid residue selection.";
				mccore::gErr(0) << std::endl;
				exit (EXIT_FAILURE);
			}
			break;
		case 's':
			{
				gstrSQLFile = optarg;
			}
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


mccore::Molecule*
loadFile (const string &filename)
{
  Molecule *molecule;
  ResidueFM rFM;
  annotate::AnnotateModelFM aFM (residueSelection, environment, &rFM);

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

bool pseudoKnots(const std::vector<annotate::Stem>& usedStems, const annotate::Stem& otherStem)
{
	bool bPseudoKnots = false;
	std::vector<annotate::Stem>::const_iterator it;
	for(it = usedStems.begin(); it != usedStems.end() && !bPseudoKnots; ++ it)
	{
		bPseudoKnots = it->pseudoKnots(otherStem);
	}
	return bPseudoKnots;
}

std::string toDotBracket(
	const annotate::AnnotationStems& aStems,
	const annotate::AnnotationChains& aChains)
{
	std::ostringstream oss;
	const std::vector< annotate::Stem >& stems = aStems.getStems();
	std::vector< annotate::Stem> usedStems;

	std::vector< annotate::Stem >::const_iterator it;
	for(it = stems.begin(); it != stems.end(); ++ it)
	{
		if(!pseudoKnots(usedStems, *it))
		{
			usedStems.push_back(*it);
		}
	}

	std::map<mccore::ResId, char> dBrackets;
	const std::list<mccore::ResId>& chain = aChains.getChain('A');
	std::list<mccore::ResId>::const_iterator itResId;
	for(itResId = chain.begin(); itResId != chain.end(); ++ itResId)
	{
		dBrackets[*itResId] = '.';
	}

	for(it = usedStems.begin(); it != usedStems.end(); ++ it)
	{
		mccore::ResId r1 = it->basePairs().front().fResId;
		mccore::ResId r2 = it->basePairs().back().fResId;
		mccore::ResId r3= it->basePairs().front().rResId;
		mccore::ResId r4 = it->basePairs().back().rResId;
		std::map<mccore::ResId, char>::iterator it1 = dBrackets.find(r1);
		std::map<mccore::ResId, char>::iterator it2 = dBrackets.find(r2);
		std::map<mccore::ResId, char>::iterator it3 = dBrackets.find(r3);
		std::map<mccore::ResId, char>::iterator it4 = dBrackets.find(r4);

		std::map<mccore::ResId, char>::iterator itWrite;
		for(itWrite = it1; itWrite != it2; itWrite ++)
		{
			itWrite->second = '(';
		}
		itWrite->second = '(';
		for(itWrite = it4; itWrite != it3; itWrite ++)
		{
			itWrite->second = ')';
		}
		itWrite->second = ')';
	}
	oss << "Number of resids : " << dBrackets.size() << std::endl;
	std::map<mccore::ResId, char>::iterator itDB;
	// Output sequence
	for(itDB = dBrackets.begin(); itDB != dBrackets.end(); ++ itDB)
	{
		oss << itDB->first << "-";
	}
	oss << std::endl;
	// Output dot-brackets
	for(itDB = dBrackets.begin(); itDB != dBrackets.end(); ++ itDB)
	{
		oss << itDB->second;
	}
	oss << std::endl;
	return oss.str();
}

int main (int argc, char *argv[])
{
	read_options (argc, argv);

	std::ofstream oSQLFile;
	annotate::AnnotatedModelSQLExporter oSQLExporter;

	if(!gstrSQLFile.empty())
	{
		oSQLFile.open (gstrSQLFile.c_str(), ios_base::out);
		if(oSQLFile.fail())
		{
			mccore::gErr (0) << PACKAGE << ": cannot open output SQL file '" << gstrSQLFile << "'." << endl;
			exit(EXIT_FAILURE);
		}
	}

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
					annotate::AnnotationInteractions annInteractions;
					annotate::AnnotationChains annChains;
					annotate::AnnotationStems annStems;
					annotate::AnnotationLinkers annLinkers;
					annotate::AnnotationLoops annLoops;
					annotate::AnnotationTertiaryPairs annTertiaryPairs;
					annotate::AnnotationTertiaryStacks annTertiaryStacks;
					annotate::AnnotationResSecondaryStructures annResSecondaryStructures;

		  			am.addAnnotation(annInteractions);
		  			am.addAnnotation(annChains);
					am.addAnnotation(annStems);
					am.addAnnotation(annLinkers);
					am.addAnnotation(annLoops);
					am.addAnnotation(annTertiaryPairs);
					am.addAnnotation(annTertiaryStacks);
					am.addAnnotation(annResSecondaryStructures);

					mccore::gOut (0) << "Annotating Model ------------------------------------------------" << endl;
					mccore::gOut (0) << filename << std::endl;
					am.annotate (gucRelationMask);
					mccore::gOut(0) << am;

					// Stems to dot bracket
					// TODO : Fix this so the chains are handled properly.
					// std::string strDotBrackets = toDotBracket(annStems, annChains);
					// mccore::gOut(0) << strDotBrackets << std::endl;

					if(oSQLFile.is_open())
					{
						oSQLExporter.addModel(am);
					}

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

	if(oSQLFile.is_open())
	{
		oSQLFile << oSQLExporter.toSQL();
		oSQLFile.close();
	}
	return EXIT_SUCCESS;
}

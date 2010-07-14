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

#include "pdb2db.h"

#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <sstream>
#include <unistd.h>

#include "mccore/Binstream.h"
#include "mccore/Exception.h"
#include "mccore/Messagestream.h"
#include "mccore/ModelFactoryMethod.h"
#include "mccore/Pdbstream.h"
#include "mccore/PropertyType.h"
#include "mccore/Relation.h"
#include "mccore/ResidueFactoryMethod.h"
#include "mccore/ResIdSet.h"
#ifdef HAVE_LIBRNAMLC__
#include "mccore/RnamlReader.h"
#endif
#include "mccore/Version.h"

// libmcannotate
#include "AlgorithmExtra.h"
#include "AnnotateModel.h"
#include "AnnotationChains.h"
#include "AnnotationInteractions.h"


#include "AnnotationStemsLoose.h"

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
const char* shortopts = "Vbe:f:hlr:vm:";
const unsigned int guiMaxLayer = 10;

std::string gNotationSymbols[] =
{
	"()ab", // 0
	"[]cd", // 1
	"<>ef", // 2
	"{}gh", // 3
	"+-ij", // 4
	"12kl", // 5
	"34mn", // 6
	"56op", // 7
	"78qr", // 8
	"90st", // 9
};

PDB2DotBracket::PDB2DotBracket(int argc, char * argv [])
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
					AnnotationStemsLoose annStems;

					am.addAnnotation(annInteractions);
					am.addAnnotation(annChains);
					am.addAnnotation(annStems);

					am.keepNucleicAcid();
					am.annotate (gucRelationMask);

					// Stems to dot bracket
					annotate::AnnotationChains::chain_map::const_iterator itChain;
					for(itChain = annChains.chains().begin(); itChain != annChains.chains().end(); ++ itChain)
					{
						mccore::gOut(0) << ">" << strPDBPrefix << ":";
						mccore::gOut(0) << itChain->first;
						mccore::gOut(0) << "|PDBID|CHAIN|SEQUENCE" << std::endl;
						mccore::gOut(0) << getSequence(am, itChain->second) << std::endl;
						std::string strDotBrackets = toDotBracketCombined(annStems, itChain->second);
						mccore::gOut(0) << strDotBrackets << std::endl;
						strDotBrackets = toDotBracketLayers(annStems, itChain->second);
						mccore::gOut(0) << strDotBrackets << std::endl;
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
}

void PDB2DotBracket::version () const
{
	mccore::Version mccorev;

	mccore::gOut(0)
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
		<< "  using " << mccorev << std::endl;
}


void PDB2DotBracket::usage () const
{
	mccore::gOut(0) << "usage: " << PACKAGE
		<< " [-bhlvV] [-e num] [-f <model number>] [-r <residue ids>] [-m <relation mask>] <structure file> ..."
		<< std::endl;
}

void PDB2DotBracket::help () const
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
		<< "  -v                be verbose" << std::endl
		<< "  -V                print the software version info" << std::endl;
}


void PDB2DotBracket::read_options (int argc, char* argv[])
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


mccore::Molecule* PDB2DotBracket::loadFile (const string &filename)
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

std::string PDB2DotBracket::getFilePrefix(const std::string& aFileName) const
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

bool PDB2DotBracket::pseudoKnots(
	const std::vector<annotate::Stem>& usedStems,
	const annotate::Stem& otherStem) const
{
	bool bPseudoKnots = false;
	std::vector<annotate::Stem>::const_iterator it;
	for(it = usedStems.begin(); it != usedStems.end() && !bPseudoKnots; ++ it)
	{
		bPseudoKnots = it->pseudoKnots(otherStem);
	}
	return bPseudoKnots;
}

std::pair<PDB2DotBracket::stem_set, PDB2DotBracket::stem_set> PDB2DotBracket::selectStems(
	const stem_set& aStems) const
{
	std::vector<annotate::Stem> stems;
	stems.insert(stems.begin(), aStems.begin(), aStems.end());
	return selectStems(stems);
}

std::pair<PDB2DotBracket::stem_set, PDB2DotBracket::stem_set> PDB2DotBracket::selectStems(
	const std::vector<annotate::Stem>& aStems) const
{
	std::pair<stem_set, stem_set> results;
	std::vector<std::set<unsigned int> > working;
	std::set<unsigned int> toTest;

	working.resize(aStems.size());
	for(unsigned int i = 0; i < aStems.size(); ++ i)
	{
		toTest.insert(i);
		for(unsigned int j = 0; j < aStems.size(); ++ j)
		{
			if(i != j)
			{
				if(aStems[i].pseudoKnots(aStems[j]) || aStems[i].overlaps(aStems[j]))
				{
					working[i].insert(j);
				}
			}
		}
	}
	std::pair<unsigned int, std::list<std::set<unsigned int> > > selection;
	selection = selectStems(toTest, working, aStems);

	for(unsigned int i = 0; i < aStems.size(); ++ i)
	{
		if(selection.second.begin()->find(i) != selection.second.begin()->end())
		{
			results.first.insert(aStems[i]);
		}else
		{
			results.second.insert(aStems[i]);

		}
	}
/*
	for(unsigned int i = 0; i < aStems.size(); ++ i)
	{
		if(0 == working[i].size())
		{
			results.first.insert(aStems[i]);
		}else
		{
			int iPseudoSize = 0;
			std::list<unsigned int>::const_iterator itCount = working[i].begin();
			for(; itCount != working[i].end(); ++ itCount)
			{
				if(results.second.end() == results.second.find(aStems[*itCount]))
				{
					iPseudoSize += aStems[*itCount].size();
				}
			}
			if(iPseudoSize < aStems[i].size())
			{
				results.first.insert(aStems[i]);
				for(itCount = working[i].begin();
					itCount != working[i].end();
					++ itCount)
				{
					results.second.insert(aStems[*itCount]);
				}
			}else if(iPseudoSize == aStems[i].size())
			{
				if(results.second.end() == results.second.find(aStems[i]))
				{
					// First encounter
					results.first.insert(aStems[i]);
					for(itCount = working[i].begin();
						itCount != working[i].end();
						++ itCount)
					{
						results.second.insert(aStems[*itCount]);
					}
				}
			}else
			{
				results.second.insert(aStems[i]);
			}
		}
	}*/
	return results;
}

std::pair<unsigned int, std::list<std::set<unsigned int> > > PDB2DotBracket::selectStems(
	const std::set<unsigned int>& aToTest,
	const std::vector<std::set<unsigned int> >& aConflicts,
	const std::vector<annotate::Stem>& aStems) const
{
	std::pair<unsigned int, std::list<std::set<unsigned int> > > returnVal;

	if(0 == aToTest.size())
	{
		std::set<unsigned int> emptySet;
		std::list<std::set<unsigned int> > emptyList;
		returnVal.first = 0;
		returnVal.second.push_back(emptySet);
	}else
	{
		unsigned int uiStem = *aToTest.begin();
		if(0 == aConflicts[uiStem].size())
		{
			// No conflicts
			std::set<unsigned int> subTest = aToTest;
			subTest.erase(uiStem);
			returnVal = selectStems(subTest, aConflicts, aStems);
			std::list<std::set<unsigned int> >::iterator it;
			for(it = returnVal.second.begin(); it != returnVal.second.end(); ++ it)
			{
				it->insert(uiStem);
			}
		}else
		{
			std::set<unsigned int> subTest1 = aToTest;
			subTest1.erase(uiStem);

			std::set<unsigned int> subTest2;
			subTest2 = annotate::SetDifference<unsigned int>(subTest1, aConflicts[uiStem]);
			std::pair<unsigned int, std::list<std::set<unsigned int> > > subVal1;
			std::pair<unsigned int, std::list<std::set<unsigned int> > > subVal2;
			subVal1 = selectStems(subTest2, aConflicts, aStems);
			subVal2 = selectStems(subTest1, aConflicts, aStems);

			// Compute the maximum value of the 2 occurrences
			unsigned int uiMaxValue = std::max(subVal1.first + aStems[uiStem].size(), subVal2.first);
			returnVal.first = uiMaxValue;
			if(uiMaxValue == (subVal1.first + aStems[uiStem].size()))
			{
				returnVal.second.insert(returnVal.second.end(), subVal1.second.begin(), subVal1.second.end());
				std::list<std::set<unsigned int> >::iterator it;
				for(it = returnVal.second.begin(); it != returnVal.second.end(); ++ it)
				{
					it->insert(uiStem);
				}
			}
			if(uiMaxValue == subVal2.first)
			{
				returnVal.second.insert(returnVal.second.end(), subVal2.second.begin(), subVal2.second.end());
			}
		}
	}
	return returnVal;
}

PDB2DotBracket::db_notation PDB2DotBracket::createDotBracket(const std::list<mccore::ResId>& aChain) const
{
	std::map<mccore::ResId, char> dBrackets;

	// Initialize the dot-bracket notation to unstructured
	std::list<mccore::ResId>::const_iterator itResId;
	for(itResId = aChain.begin(); itResId != aChain.end(); ++ itResId)
	{
		dBrackets[*itResId] = '.';
	}
	return dBrackets;
}

std::string PDB2DotBracket::toDotBracketLayers(
	const AnnotationStemsLoose& aStems,
	const std::list<mccore::ResId>& aChain) const
{
	const char cChain = aChain.begin()->getChainId();
	std::list<std::set<annotate::Stem> > layers = splitLayers(cChain, aStems);
	return toDotBracket(aChain, layers);
}

std::string PDB2DotBracket::toDotBracketCombined(
	const AnnotationStemsLoose& aStems,
	const std::list<mccore::ResId>& aChain) const
{
	const char cChain = aChain.begin()->getChainId();
	std::list<std::set<annotate::Stem> > layers = splitLayers(cChain, aStems);
	std::ostringstream oss;
	db_notation dBrackets = createDotBracket(aChain);

	unsigned int uiLayer = 0;
	std::list<std::set<annotate::Stem> >::const_iterator it;
	for(it = layers.begin(); it != layers.end() && uiLayer < guiMaxLayer; ++ it)
	{
		applyStems(dBrackets, *it, uiLayer);
		uiLayer ++;
	}

	// Output dot-brackets
	std::map<mccore::ResId, char>::iterator itDB;
	for(itDB = dBrackets.begin(); itDB != dBrackets.end(); ++ itDB)
	{
		oss << itDB->second;
	}
	return oss.str();
}

std::string PDB2DotBracket::toDotBracket(
	const std::list<mccore::ResId>& aChain,
	const std::list<std::set<annotate::Stem> >& aLayers) const
{
	std::ostringstream oss;

	unsigned int uiLayer = 0;
	std::list<std::set<annotate::Stem> >::const_iterator it;
	for(it = aLayers.begin(); it != aLayers.end() && uiLayer < guiMaxLayer; ++ it)
	{
		std::list<mccore::ResId>::const_iterator itResId;
		db_notation dBrackets = createDotBracket(aChain);

		applyStems(dBrackets, *it, uiLayer);

		// Output dot-brackets
		if(it != aLayers.begin())
		{
			oss << std::endl;
		}
		std::map<mccore::ResId, char>::iterator itDB;
		for(itDB = dBrackets.begin(); itDB != dBrackets.end(); ++ itDB)
		{
			oss << itDB->second;
		}
		uiLayer ++;
	}
	return oss.str();
}

std::list<PDB2DotBracket::stem_set> PDB2DotBracket::splitLayers(
	const char acChain,
	const AnnotationStemsLoose& aStems) const
{
	std::list<stem_set> layers;
	const std::vector<annotate::Stem>& stems = aStems.getStems();
	std::set<annotate::Stem> usedStems;
	std::vector< annotate::Stem>::const_iterator it;
	for(it = stems.begin(); it != stems.end(); ++ it)
	{
		char cStemChain = it->basePairs().front().fResId.getChainId();
		if(acChain == cStemChain)
		{
			usedStems.insert(*it);
		}
	}

	// Select the non pseudo-knotted stems
	std::pair<stem_set, stem_set> selectedStems;
	selectedStems.second = usedStems;

	while(0 < selectedStems.second.size())
	{
		selectedStems = selectStems(selectedStems.second);
		layers.push_back(selectedStems.first);
	}
	return layers;
}

void PDB2DotBracket::applyStems(
	PDB2DotBracket::db_notation& aDBNotation,
	const std::set< annotate::Stem>& aStems,
	unsigned int auiLevel) const
{
	assert(auiLevel < guiMaxLayer);
	std::set< annotate::Stem >::const_iterator it;
	for(it = aStems.begin(); it != aStems.end(); ++ it)
	{
		const std::vector< annotate::BasePair>& pairs = it->basePairs();
		char cOpen = gNotationSymbols[auiLevel][0];
		char cClose = gNotationSymbols[auiLevel][1];
		if(it->getOrientation() == annotate::Stem::ePARALLEL)
		{
			cOpen = gNotationSymbols[auiLevel][2];
			cClose = gNotationSymbols[auiLevel][3];
		}

		std::vector< annotate::BasePair>::const_iterator itPair;
		for(itPair = pairs.begin(); itPair != pairs.end(); ++ itPair)
		{
			applyPair(aDBNotation, *itPair, cOpen, cClose);
		}
	}
}

void PDB2DotBracket::applyPair(
	PDB2DotBracket::db_notation& aDBNotation,
	const annotate::BasePair& aPair,
	const char& acOpen,
	const char& acClose) const
{
	std::map<mccore::ResId, char>::iterator it1 = aDBNotation.find(aPair.fResId);
	std::map<mccore::ResId, char>::iterator it2 = aDBNotation.find(aPair.rResId);

	if(it1->second == '.' && it2->second == '.')
	{
		it1->second = acOpen;
		it2->second = acClose;
	}
}

std::string PDB2DotBracket::getSequence(
	const annotate::AnnotateModel& am,
	const std::list<mccore::ResId>& aChain) const
{
	std::ostringstream oss;
	std::list<mccore::ResId>::const_iterator it = aChain.begin();
	for(; it != aChain.end(); ++ it)
	{
		annotate::AnnotateModel::const_iterator itRes = am.find(*it);
		assert(itRes != am.end());
		oss << mccore::Pdbstream::stringifyResidueType (itRes->getType ());
	}
	return oss.str();
}

int main(int argc, char* argv[])
{
	PDB2DotBracket theApp(argc, argv);
	return EXIT_SUCCESS;
}

//                              -*- Mode: C++ -*-
// pdb2db.cc
// Copyright © 2001-10 Laboratoire d'ingénierie des ARN.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Thu Jul 8 10:00:00 2010


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
ResIdSet residueSelection;
const char* shortopts = "Vbe:f:hlr:vm:";
const unsigned int guiMaxLayer = 10;
unsigned int guiMaxPerfectSearch = 24;

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
		mccore::Molecule *molecule;
		mccore::Molecule::iterator molIt;
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
					unsigned char ucRelationMask =
						mccore::Relation::adjacent_mask
						| mccore::Relation::pairing_mask;
					am.annotate(ucRelationMask);

					// Stems to dot bracket
					annotate::AnnotationChains::chain_map::const_iterator itChain;
					for(itChain = annChains.chains().begin(); itChain != annChains.chains().end(); ++ itChain)
					{
						std::cout << ">" << strPDBPrefix << ":";
						std::cout << itChain->first;
						std::cout << "|PDBID|CHAIN|SEQUENCE" << std::endl;
						std::cout << getSequence(am, itChain->second) << std::endl;
						std::string strDotBrackets;
						strDotBrackets = toDotBracketCombined(annStems, itChain->second);
						std::cout << strDotBrackets << std::endl;
						strDotBrackets = toDotBracketLayers(annStems, itChain->second);
						std::cout << strDotBrackets << std::endl;
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
		<< " [-hlvV] [-f <model number>] [-m <relation mask>] <structure file> ..."
		<< std::endl;
}

void PDB2DotBracket::help () const
{
	mccore::gOut (0)
		<< "This program annotate structures (and more)." << std::endl
		<< "  -s 			    " << std::endl
		<< "  -f model number   model to annotate (default : all models) " << std::endl
		<< "  -h                print this help" << std::endl
		<< "  -l                be more verbose (log)" << std::endl
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

	if(aStems.size() < guiMaxPerfectSearch)
	{
		std::vector<std::set<unsigned int> > conflicts = computeConflicts(aStems);

		// Set of the stem indices that needs testing
		std::set<unsigned int> toTest;
		for(unsigned int i = 0; i < aStems.size(); ++ i)
		{
			toTest.insert(i);
		}

		std::pair<unsigned int, std::list<std::set<unsigned int> > > selection;
		selection = selectStemsRecursive(toTest, conflicts, aStems);

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
	}else
	{
		std::set<annotate::Stem> stems;
		stems.insert(aStems.begin(), aStems.end());
		results = splitImbrication(stems);
	}
	return results;
}

std::vector<std::set<unsigned int> > PDB2DotBracket::computeConflicts(
	const std::vector<annotate::Stem>& aStems) const
{
	std::vector<std::set<unsigned int> > conflicts;
	conflicts.resize(aStems.size());
	for(unsigned int i = 0; i < aStems.size(); ++ i)
	{
		for(unsigned int j = 0; j < aStems.size(); ++ j)
		{
			if(i != j)
			{
				if(aStems[i].pseudoKnots(aStems[j]) || aStems[i].overlaps(aStems[j]))
				{
					conflicts[i].insert(j);
				}
			}
		}
	}
	return conflicts;
}

std::pair<PDB2DotBracket::stem_set, PDB2DotBracket::stem_set> PDB2DotBracket::splitConflicts(
	const PDB2DotBracket::stem_set& aStems) const
{
	std::pair<PDB2DotBracket::stem_set, PDB2DotBracket::stem_set> returnVal;
	std::vector<std::set<unsigned int> > conflicts;
	std::vector<annotate::Stem> stems;
	stems.insert(stems.begin(), aStems.begin(), aStems.end());

	conflicts = computeConflicts(stems);
	for(unsigned int i = 0; i < stems.size(); ++ i)
	{
		if(0 < conflicts[i].size())
		{
			returnVal.second.insert(stems[i]);
		}else
		{
			returnVal.first.insert(stems[i]);
		}
	}
	return returnVal;
}

std::pair<unsigned int, std::list<std::set<unsigned int> > > PDB2DotBracket::selectStemsRecursive(
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
			returnVal = selectStemsRecursive(subTest, aConflicts, aStems);
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
			subVal1 = selectStemsRecursive(subTest2, aConflicts, aStems);
			subVal2 = selectStemsRecursive(subTest1, aConflicts, aStems);

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

unsigned int PDB2DotBracket::findMaxConflictStem(
	const std::set<unsigned int>& aToTest,
	const std::vector<std::set<unsigned int> >& aConflicts,
	const std::vector<annotate::Stem>& aStems) const
{
	unsigned int uiStem = *aToTest.begin();
	unsigned int uiMaxConflicts = aConflicts[uiStem].size();
	std::set<unsigned int>::const_iterator it;
	for(it = aToTest.begin(); it != aToTest.end(); ++ it)
	{
		if(aConflicts[*it].size() > uiMaxConflicts)
		{
			uiMaxConflicts = aConflicts[*it].size();
			uiStem = *it;
		}
	}
	return uiStem;
}

std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > PDB2DotBracket::splitImbrication(
	const std::set<annotate::Stem>& aStems) const
{
	std::set<annotate::Stem> work = aStems;
	std::set<annotate::Stem> toSelect;
	std::set<annotate::Stem> toReserve;

	while(!work.empty())
	{
		if(1 == work.size())
		{
			// Only one left, keep it
			toSelect.insert(*work.begin());
			work.clear();
		}else
		{
			std::vector<annotate::Stem> stemVector;
			stemVector.insert(stemVector.begin(), work.begin(), work.end());
			std::vector<std::set<unsigned int> > conflicts = computeConflicts(stemVector);
			unsigned int uiMaxStem = 0;
			unsigned int uiMaxCount = 0;
			for(unsigned int i = 0; i < stemVector.size(); ++ i)
			{
				unsigned int uiCount = 0;
				std::set<unsigned int>::const_iterator it;
				for(it = conflicts[i].begin(); it != conflicts[i].end(); ++ it)
				{
					uiCount += stemVector[*it].size();
					if(uiMaxCount < uiCount)
					{
						uiMaxCount = uiCount;
						uiMaxStem = i;
					}
				}
			}
			// Remove max stem
			toReserve.insert(stemVector[uiMaxStem]);
			work.erase(stemVector[uiMaxStem]);
			for(unsigned int j = 0; j < conflicts.size(); ++ j)
			{
				conflicts[j].erase(uiMaxStem);
				if(0 == conflicts[j].size())
				{
					work.erase(stemVector[j]);
					toSelect.insert(stemVector[j]);
				}
			}
		}
	}
	return std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> >(toSelect, toReserve);
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
	std::ostringstream oss;
	std::list<stem_set> layers;

	// Select the non pseudo-knotted stems
	std::pair<stem_set, stem_set> selectedStems;
	selectedStems.second = keepChain(cChain, aStems.getStems());
	stem_set assignedStems;

	while(0 < selectedStems.second.size())
	{
		// Split the conflicting stems from the simple ones
		std::pair<stem_set, stem_set> conflictSplit = splitConflicts(selectedStems.second);
		selectedStems.second = conflictSplit.second;

		// Apply selection on conflicting stems
		selectedStems = selectStems(selectedStems.second);

		// Keep the selected stems and the unconflicted ones for this layer
		selectedStems.first.insert(conflictSplit.first.begin(), conflictSplit.first.end());
		assignedStems.insert(selectedStems.first.begin(), selectedStems.first.end());
		layers.push_back(selectedStems.first);

		// Remove the already assigned residues from the remaining stems
		unsigned int i = 0;
		stem_set::const_iterator it = selectedStems.second.begin();
		for(;it != selectedStems.second.end(); ++ it)
		{
			std::set<annotate::Stem> test;
			test.insert(*it);
			test = cutStems(test, assignedStems);
			++ i;
		}
		selectedStems.second = cutStems(selectedStems.second, assignedStems);
	}

	db_notation dBrackets = createDotBracket(aChain);
	unsigned int uiLayer = 0;
	std::list<stem_set>::const_iterator it;
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

std::string PDB2DotBracket::toDotBracket(
	const std::list<mccore::ResId>& aChain,
	const stem_set& aLayer) const
{
	std::ostringstream oss;

	db_notation dBrackets = createDotBracket(aChain);
	applyStems(dBrackets, aLayer, 0);

	std::map<mccore::ResId, char>::iterator itDB;
	for(itDB = dBrackets.begin(); itDB != dBrackets.end(); ++ itDB)
	{
		oss << itDB->second;
	}
	return oss.str();
}

std::list<PDB2DotBracket::stem_set> PDB2DotBracket::splitLayers(
	const char acChain,
	const AnnotationStemsLoose& aStems) const
{
	std::list<stem_set> layers;
	std::set<annotate::Stem> usedStems;
	usedStems = keepChain(acChain, aStems.getStems());

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

PDB2DotBracket::stem_set PDB2DotBracket::keepChain(
	const char acChain,
	const std::vector<annotate::Stem>& aStems) const
{
	stem_set keptStems;
	std::set<annotate::Stem> usedStems;
	usedStems.insert(aStems.begin(), aStems.end());
	return keepChain(acChain, usedStems);
}

PDB2DotBracket::stem_set PDB2DotBracket::keepChain(
	const char acChain,
	const std::set<annotate::Stem>& aStems) const
{
	stem_set keptStems;
	std::set<annotate::Stem>::const_iterator it;
	for(it = aStems.begin(); it != aStems.end(); ++ it)
	{
		char cStemChain = it->basePairs().front().fResId.getChainId();
		if(acChain == cStemChain)
		{
			keptStems.insert(*it);
		}
	}
	return keptStems;
}

void PDB2DotBracket::applyStems(
	PDB2DotBracket::db_notation& aDBNotation,
	const std::set<annotate::Stem>& aStems,
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

PDB2DotBracket::stem_set PDB2DotBracket::cutStem(
	const annotate::Stem& aStem,
	const stem_set& aStems) const
{
	std::set<annotate::Stem> stems;
	std::set<annotate::Stem>::const_iterator itStem;

	std::set<annotate::BasePair> stemPairs;
	stemPairs.insert(aStem.basePairs().begin(), aStem.basePairs().end());

	std::set<mccore::ResId> removedResId;
	for(itStem = aStems.begin(); itStem != aStems.end(); ++ itStem)
	{
		std::vector<annotate::BasePair>::const_iterator it;
		for(it = itStem->basePairs().begin(); it != itStem->basePairs().end(); ++ it)
		{
			removedResId.insert(it->fResId);
			removedResId.insert(it->rResId);
		}
	}

	std::set<annotate::BasePair> remainingPairs;
	std::set<annotate::BasePair>::const_iterator itPair;
	for(itPair = stemPairs.begin(); itPair != stemPairs.end(); ++ itPair)
	{
		if((removedResId.end() == removedResId.find(itPair->fResId)) && (removedResId.end() == removedResId.find(itPair->rResId)))
		{
			remainingPairs.insert(*itPair);
		}
	}

	annotate::Stem stem;
	for(itPair = remainingPairs.begin(); itPair != remainingPairs.end(); ++ itPair)
	{
		if(!stem.continues(*itPair))
		{
			stems.insert(stem);
			stem.clear();
			stem.push_back(*itPair);
		}else
		{
			stem.push_back(*itPair);
		}
	}
	if(0 < stem.size())
	{
		stems.insert(stem);
		stem.clear();
	}

	return stems;
}

PDB2DotBracket::stem_set PDB2DotBracket::cutStems(
	const stem_set& aToCutStems,
	const stem_set& aStems) const
{
	stem_set stems;
	stem_set cut;
	stem_set::const_iterator it;
	for(it = aToCutStems.begin(); it != aToCutStems.end(); ++ it)
	{
		annotate::Stem toCutStem = *it;
		cut = cutStem(toCutStem, aStems);
		stems.insert(cut.begin(), cut.end());
	}
	return stems;
}

int main(int argc, char* argv[])
{
	PDB2DotBracket theApp(argc, argv);
	return EXIT_SUCCESS;
}

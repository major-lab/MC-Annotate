//                              -*- Mode: C++ -*-
// mcssc.cc
// Copyright © 2009 Laboratoire de Biologie Informatique et Théorique.
//                  Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Aug 18 11:38:07 2009
// $Revision: 1 $
// $Id: mcssc.cc 1 2009-09-18 11:38:07 blanchmf $



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cerrno>
#include <cstdlib>
#include <string>
#include <sstream>
#include <unistd.h>
#include <cassert>

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
#include "AnnotationStems.h"
#include "BaseLink.h"
#include "BasePair.h"
#include "Cycle.h"
#include "StringTable.h"

bool binary = false;
bool oneModel = false;
unsigned int guiModelNumber = 0;  // 1 based vector identifier, 0 means all
unsigned char gucRelationMask =
	mccore::Relation::adjacent_mask
	| mccore::Relation::pairing_mask
	| mccore::Relation::stacking_mask
	| mccore::Relation::backbone_mask;
const char* shortopts = "Vbf:hlvm:";

// PROTOTYPES ------------------------------------------------------------------
std::string getResIdsString(const annotate::Cycle& aCycle);

void version ()
{
	mccore::Version mccorev;

	mccore::gOut (0) << PACKAGE << " " << VERSION << " (" << __DATE__ << ")" << std::endl
	   << "  using " << mccorev << std::endl;
}


void usage ()
{
	mccore::gOut (0) << "usage: " << PACKAGE
	   << " [-bhlvV] [-f <model number>] [-m <relation mask>] <structure file> ..."
	   << std::endl;
}

void help ()
{
	mccore::gOut(0)
		<< "This program annotate secondary structures and produces output the cycles found." << endl
		<< "  -b                read binary files instead of pdb files" << endl
		<< "  -m <mask>         annotation mask string: any combination of 'A' (adjacent), 'S' (stacking), 'P' (pairing) and 'B' (backbone). (default: all)" << endl
		<< "  -f model number   model to print" << endl
		<< "  -h                print this help" << endl
		<< "  -l                be more verbose (log)" << endl
		<< "  -v                be verbose" << endl
		<< "  -V                print the software version info" << endl;
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
	annotate::AnnotateModelFM aFM (residueSelection, 0, &rFM);

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

void cleanString(std::string& aString, const char& aChar)
{
	std::string::iterator it = aString.begin();
	while(it != aString.end())
	{
		if(*it == aChar)
		{
			it = aString.erase(it);
		}
		else
		{
			++ it;
		}
	}
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
	cleanString(filename, ' ');
	return filename;
}

std::set<annotate::Cycle> computeStemCycles(
	const mccore::GraphModel& aModel,
	const annotate::Stem& aStem)
{
	std::set<annotate::Cycle> cycles;

	std::vector< annotate::BasePair >::const_iterator it;
	std::vector< annotate::BasePair >::const_iterator itPrev;
	for(it = aStem.basePairs().begin();	it != aStem.basePairs().end();	++ it)
	{
		if(it != aStem.basePairs().begin())
		{
			annotate::Cycle::interactions_set interactions;

			annotate::BasePair bp1(
				itPrev->first, itPrev->fResId,
				itPrev->second, itPrev->rResId);
			annotate::BasePair bp2(
				it->first, it->fResId,
				it->second, it->rResId);
			annotate::BaseLink bl1(
				itPrev->first, itPrev->fResId,
				it->first, it->fResId);
			annotate::BaseLink bl2(
				it->second, it->rResId,
				itPrev->second, itPrev->rResId);

			interactions.insert(bp1);
			interactions.insert(bp2);
			interactions.insert(bl1);
			interactions.insert(bl2);

			// Insert the residues in the model
			annotate::Cycle cycle(interactions);

			// clearInteractionsSet(interactions);
			cycles.insert(cycle);
		}
		itPrev = it;
	}
	return cycles;
}

annotate::Cycle computeLoopCycle(
	const annotate::AnnotateModel &aModel,
	const annotate::Loop& aLoop)
{
	assert(0 < aLoop.linkers().size());

	annotate::Cycle::interactions_set interactions;

	interactions = aLoop.getBaseInteractions();

	// Insert the residues in the model
	assert(0 < interactions.size());
	annotate::Cycle cycle(interactions);

	return cycle;
}

std::string debugCycleString(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;

	std::set<annotate::BaseInteraction> inters = aCycle.getBaseInteractions();

	for(std::set<annotate::BaseInteraction>::const_iterator it = inters.begin(); it != inters.end(); ++ it)
	{
		if(it != inters.begin())
		{
			oss << ",";
		}
		oss << "{" << it->fResId << "-" << it->rResId << "}";
	}
	return oss.str();
}

std::set<annotate::BasePair> getDividingPairs(
	const annotate::Cycle& aCycle,
	const std::set<annotate::BasePair>& aPairs)
{
	std::set<annotate::BaseInteraction> cycleInteractions = aCycle.getBaseInteractions();
	std::set<annotate::BasePair> dividingPairs;
	std::list<mccore::ResId> resIds = aCycle.resIds();
	std::set<annotate::BasePair>::const_iterator itPair;
	for(itPair = aPairs.begin(); itPair != aPairs.end(); ++ itPair)
	{
		if(resIds.end() != std::find(resIds.begin(), resIds.end(), itPair->fResId)
			&& resIds.end() != std::find(resIds.begin(), resIds.end(), itPair->rResId))
		{
			annotate::BaseInteraction search(0, itPair->fResId, 0, itPair->rResId);
			if(cycleInteractions.end() == cycleInteractions.find(search))
			{
				dividingPairs.insert(*itPair);
			}
		}
	}

	return dividingPairs;
}

bool isValidNCM(const annotate::Cycle& aCycle)
{
	bool bIsValid = false;

	switch(aCycle.getType())
	{
	case annotate::Cycle::eLOOSE:
		std::cout << "eLOOSE" << std::endl;
		bIsValid = true;
		break;
	case annotate::Cycle::eLOOP:
		std::cout << "eLOOP" << std::endl;
		bIsValid = true;
		break;
	case annotate::Cycle::e2STRANDS_TRIANGLE:
		// Don't support the triangles
		std::cout << "e2STRANDS_TRIANGLE" << std::endl;
		break;
	case annotate::Cycle::e2STRANDS_PARALLEL:
		std::cout << "e2STRANDS_PARALLEL" << std::endl;
		bIsValid = true;
		break;
	case annotate::Cycle::e2STRANDS_ANTIPARALLEL:
		std::cout << "e2STRANDS_ANTIPARALLEL" << std::endl;
		bIsValid = true;
		break;
	case annotate::Cycle::eMULTIBRANCH:
		std::cout << "eMULTIBRANCH" << std::endl;
		bIsValid = true;
		break;
	}
	return bIsValid;
}

bool isNewKnowledge(const annotate::Cycle& aCycle)
{
	bool bIsNew = false;

	std::cout << debugCycleString(aCycle) << std::endl;
	std::cout << "Ordered residues : ";
	std::list<mccore::ResId>::const_iterator it = aCycle.resIds().begin();
	for(; it != aCycle.resIds().end(); ++ it)
	{
		if(it != aCycle.resIds().begin())
		{
			std::cout << ",";
		}
		std::cout << *it;
	}
	std::cout << std::endl;

	assert(isValidNCM(aCycle));

	switch(aCycle.getType())
	{
	case annotate::Cycle::eLOOSE:
		bIsNew = true;
		break;
	case annotate::Cycle::eLOOP:
		bIsNew = (6 < aCycle.resIds().size());
		break;
	case annotate::Cycle::e2STRANDS_PARALLEL:
		bIsNew = (8 < aCycle.resIds().size());
		break;
	case annotate::Cycle::e2STRANDS_ANTIPARALLEL:
		bIsNew = (8 < aCycle.resIds().size());
		break;
	case annotate::Cycle::eMULTIBRANCH:
		bIsNew = true;
		break;
	default:
		bIsNew = false;
	}
	return bIsNew;
}

annotate::Cycle makeCycle(
	const annotate::Cycle& aCycle,
	const std::list<mccore::ResId>& aResidues,
	const annotate::BasePair& aPair)
{
	std::set<mccore::ResId> residues;
	residues.insert(aResidues.begin(), aResidues.end());
	annotate::Cycle::interactions_set cycleInteractions;
	annotate::Cycle::interactions_set interactions = aCycle.getInteractions();
	annotate::Cycle::interactions_set::const_iterator itInter;
	for(itInter = interactions.begin(); itInter != interactions.end(); ++ itInter)
	{
		assert(itInter->type() != annotate::BaseInteraction::eUNKNOWN);
		if( residues.end() != residues.find(itInter->fResId)
			&& residues.end() != residues.find(itInter->rResId))
		{
			cycleInteractions.insert(*itInter);
		}
	}
	cycleInteractions.insert(aPair);

	annotate::Cycle cycle(cycleInteractions);

	return cycle;
}

std::list<annotate::Cycle> divideCycle(
	const annotate::Cycle &aCycle,
	const annotate::BasePair& aPair)
{
	std::list<annotate::Cycle> dividedCycles;
	std::list<mccore::ResId> orderedResIds = aCycle.resIds();
	std::list<mccore::ResId> cycle1ResIds;
	std::list<mccore::ResId> cycle2ResIds;
	std::list<mccore::ResId>::const_iterator itResId;

	// Before first cutting point
	for(itResId = orderedResIds.begin(); itResId != orderedResIds.end(); ++ itResId)
	{
		cycle1ResIds.push_back(*itResId);
		if(*itResId == aPair.fResId || *itResId == aPair.rResId)
		{
			cycle2ResIds.push_back(*itResId);
			++ itResId;
			break;
		}
	}

	// Between cutting points
	for(; itResId != orderedResIds.end(); ++ itResId)
	{
		cycle2ResIds.push_back(*itResId);
		if(*itResId == aPair.fResId || *itResId == aPair.rResId)
		{
			cycle1ResIds.push_back(*itResId);
			++ itResId;
			break;
		}
	}

	// After cutting points
	for(; itResId != orderedResIds.end(); ++ itResId)
	{
		cycle1ResIds.push_back(*itResId);
	}

	// Create cycles
	annotate::Cycle cycle1 = makeCycle(aCycle, cycle1ResIds, aPair);
	annotate::Cycle cycle2 = makeCycle(aCycle, cycle2ResIds, aPair);

	// Verify that they are valid NCM cycle
	if(isValidNCM(cycle1))
	{
		dividedCycles.push_back(cycle1);
	}
	if(isValidNCM(cycle2))
	{
		dividedCycles.push_back(cycle2);
	}

	// If there is no valid cycle, simply return the one passed in parameter
	if(0 == dividedCycles.size())
	{
		dividedCycles.push_back(aCycle);
	}

	return dividedCycles;
};

std::list<annotate::Cycle> divideCycle(
	const annotate::Cycle &aCycle,
	const std::set<annotate::BasePair>& aPairs)
{
	std::list<annotate::Cycle> dividedCycles;
	std::set<annotate::BasePair> pairs = aPairs;
	while(!pairs.empty())
	{
		annotate::BasePair pair = *pairs.begin();
		pairs.erase(pairs.begin());
		std::list<annotate::Cycle> divided = divideCycle(aCycle, pair);
		if(1 < divided.size())
		{
			std::set<annotate::BasePair> dividingPairs1 = getDividingPairs(divided.front(), pairs);
			std::list<annotate::Cycle> div1 = divideCycle(divided.front(), dividingPairs1);

			std::set<annotate::BasePair> dividingPairs2 = getDividingPairs(divided.back(), pairs);
			std::list<annotate::Cycle> div2 = divideCycle(divided.back(), dividingPairs2);

			dividedCycles.insert(dividedCycles.end(), div1.begin(), div1.end());
			dividedCycles.insert(dividedCycles.end(), div2.begin(), div2.end());
			break;
		}
	}

	// If we didn't manage to divide anything, return the cycle in parameter
	if(dividedCycles.empty())
	{
		dividedCycles.push_back(aCycle);
	}

	return dividedCycles;
}

std::set<annotate::Cycle> computeSecondaryStructureCycles(
	const annotate::AnnotateModel &aModel)
{
	std::set<annotate::Cycle> cycles;
	const annotate::AnnotationInteractions* pInteractions = NULL;
	const annotate::AnnotationStems* pStems = NULL;
	const annotate::AnnotationLoops* pLoops = NULL;

	pStems = aModel.getAnnotation<annotate::AnnotationStems>();
	pLoops = aModel.getAnnotation<annotate::AnnotationLoops>();
	pInteractions = aModel.getAnnotation<annotate::AnnotationInteractions>();

	std::vector<annotate::Stem>::const_iterator itStem;
	for(itStem = pStems->getStems().begin();
		itStem != pStems->getStems().end();
		++ itStem)
	{
		std::set<annotate::Cycle> stemCycles = computeStemCycles(aModel, *itStem);
		cycles.insert(stemCycles.begin(), stemCycles.end());
	}

	std::vector<annotate::Loop>::const_iterator itLoop;
	for(itLoop = pLoops->getLoops().begin();
		itLoop != pLoops->getLoops().end();
		++ itLoop)
	{
		annotate::Cycle cycle = computeLoopCycle(aModel, *itLoop);
		std::set<annotate::BasePair> pairs;
		std::set<annotate::BasePair> interPairs;
		interPairs.insert(pInteractions->pairs().begin(), pInteractions->pairs().end());
		pairs = getDividingPairs(cycle, interPairs);
		std::list<annotate::Cycle> dividedCycles = divideCycle(cycle, pairs);
		cycles.insert(dividedCycles.begin(), dividedCycles.end());
		cycles.insert(cycle);
	}
	return cycles;
}

std::string cycleProfileStrandString(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;

	// Describe the profile
	std::vector<unsigned int>::const_iterator it = aCycle.profile().begin();
	for(; it != aCycle.profile().end(); ++ it)
	{
		if(it != aCycle.profile().begin())
		{
			oss << "_";
		}
		oss << *it;
	}
	return oss.str();
}

std::string getModelString(const unsigned int auiModel)
{
	std::ostringstream oss;
	oss << auiModel;
	return oss.str();
}

std::string getResIdsString(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;

	std::list<mccore::ResId>::const_iterator it;
	for(it = aCycle.resIds().begin(); it != aCycle.resIds().end(); ++ it)
	{
		if(it != aCycle.resIds().begin())
		{
			oss << "-";
		}
		oss << *it;
	}
	return oss.str();
}

std::string getResiduesString(
	const annotate::AnnotateModel &aModel,
	const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	std::list<mccore::ResId>::const_iterator it;
	for(it = aCycle.resIds().begin(); it != aCycle.resIds().end(); ++ it)
	{
		if(it != aCycle.resIds().begin())
		{
			oss << "-";
		}
		GraphModel::const_iterator itRes = aModel.find(*it);
		assert(itRes != aModel.end());
		oss << mccore::Pdbstream::stringifyResidueType(itRes->getType());
	}

	return oss.str();
}

std::string cycleProfileString(const annotate::Cycle& aCycle)
{
	std::ostringstream oss;
	// Describe the profile
	annotate::Cycle::enType eCycleType = aCycle.getType();
	switch(eCycleType)
	{
		case annotate::Cycle::eLOOSE:
			oss << aCycle.resIds().size();
			oss << "_L";
			break;
		case annotate::Cycle::e2STRANDS_PARALLEL:
			oss << cycleProfileStrandString(aCycle) << "_p";
			break;
		case annotate::Cycle::eMULTIBRANCH:
			oss << aCycle.resIds().size();
			oss << "_M";
			break;
		default:
			oss << cycleProfileStrandString(aCycle);
			break;
	}

	return oss.str();
}

int main (int argc, char *argv[])
{
	annotate::StringTable stringTable(6);
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
					std::cout << filename << std::endl;
					annotate::AnnotateModel &am = (annotate::AnnotateModel&) *molIt;
					am.name(getPdbFileName(filename));
					annotate::AnnotationInteractions annInteractions;
					annotate::AnnotationChains annChains;
					annotate::AnnotationStems annStems;
					annotate::AnnotationLinkers annLinkers;
					annotate::AnnotationLoops annLoops;

		  			am.addAnnotation(annInteractions);
		  			am.addAnnotation(annChains);
					am.addAnnotation(annStems);
					am.addAnnotation(annLinkers);
					am.addAnnotation(annLoops);

					am.annotate (gucRelationMask);

					std::set<annotate::Cycle> cycles;
					cycles = computeSecondaryStructureCycles(am);
					std::set<annotate::Cycle>::const_iterator itCycle;
					for( itCycle = cycles.begin(); itCycle != cycles.end(); ++ itCycle)
					{
						if(isNewKnowledge(*itCycle))
						{
							std::vector<string>& tableRow = stringTable.addRow();
							tableRow[0] = getPdbFileName(filename);
							tableRow[1] = getModelString(uiModel);
							tableRow[2] = cycleProfileString(*itCycle);
							tableRow[3] = cycleProfileString(*itCycle);
							tableRow[4] = getResIdsString(*itCycle);
							tableRow[5] = getResiduesString(am, *itCycle);
						}
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
	mccore::gOut(0) << stringTable.toString(" : ");

	return EXIT_SUCCESS;
}

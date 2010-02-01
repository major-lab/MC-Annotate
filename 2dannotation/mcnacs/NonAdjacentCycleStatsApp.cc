/*
 * NonAdjacentCycleStatsApp.cc
 *
 *  Created on: Dec 15, 2009
 *      Author: blanchmf
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "NonAdjacentCycleStatsApp.h"

#include "Interaction.h"
#include "StringUtil.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

const char* shortopts = "Vhup:c:d:";

NonAdjacentCycleStatsApp::NonAdjacentCycleStatsApp(int argc, char * argv [])
{
	mbSplitInteractions = false;
	mbRemoveComposites = false;

	readOptions(argc, argv);

	readInteractionsFile();
	readCyclesFile();

	// Identify the non adjacent cycles ( cycles containing non-adjacent interactions )
	mNACycles = getNonAdjacentCycles();

	// Identify the cycles of the secondary structure
	mSecondaryStructureCycles = annotate::SetDifference(mCycles, mNACycles);

	// Identify the connections of the cycle
	// TODO: Decide if this needs to be removed.
	mConnectedCycles = this->getConnectedCycles(mNACycles);
}

void NonAdjacentCycleStatsApp::version () const
{
	std::cout
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")"
		<< std::endl;
}

void NonAdjacentCycleStatsApp::usage () const
{
	std::cout
	 	<< "usage: " << PACKAGE
		<< " [-fhuvV] [-d <output directory>] -p <distant pairs file> -c <cycles file>"
		<< std::endl;
}

void NonAdjacentCycleStatsApp::help () const
{
	std::cout
		<< "This program read cycles and 3D interactions and compute stats on them."
		<< std::endl
		<< "  -d	specify an output directory for additional annotation" << std::endl
		<< "  -f	filter out composite cycles" << std::endl
		<< "  -h	print this help" << std::endl
		<< "  -u	split the non-adjacent cycles into sets of unique interacting cycle per strand" << std::endl
		<< "  -p	file containing the non-adjacent interacting pairs" << std::endl
		<< "  -c	file containing the identified cycles" << std::endl
		<< "  -V	print the software version info" << std::endl;
}

void NonAdjacentCycleStatsApp::readOptions (int argc, char* argv[])
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
		case 'c':
		{
			mstrCyclesFile = optarg;
			break;
		}
		case 'f':
		{
			mbRemoveComposites = true;
			break;
		}
		case 'p':
		{
			mstrInteractionsFile = optarg;
			break;
		}
		case 'd':
		{
			mstrOutputDirectory = optarg;
		}
		case 'u':
		{
			mbSplitInteractions = true;
			break;
		}
		case 'h':
			usage ();
			help ();
			exit (EXIT_SUCCESS);
			break;
		default:
			usage ();
			exit (EXIT_FAILURE);
		}
	}

	if(mstrInteractionsFile.empty())
	{
		usage ();
		exit (EXIT_FAILURE);
	} else if(mstrCyclesFile.empty())
	{
		usage ();
		exit (EXIT_FAILURE);
	}
}

/**
 * readInteractionsFile
 * @brief Read interactions from a file.
 */
void NonAdjacentCycleStatsApp::readInteractionsFile()
{
	std::ifstream infile;
	mInteractions.clear();

	infile.open (mstrInteractionsFile.c_str(), std::ifstream::in);
	std::string strLine;

	while (std::getline(infile, strLine).good())
	{
		annotate::cleanString(strLine, ' ');

		std::size_t sep = strLine.rfind('-');
		std::string strRes2 = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);

		sep = strLine.rfind(':');
		std::string strRes1 = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);

		sep = strLine.rfind(':');
		std::string strModel = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);
		unsigned int uiModel = atol(strModel.c_str());

		std::string strPDBFile = strLine;

		annotate::InteractionInfo info(strPDBFile, uiModel, strRes1, strRes2);
		mInteractions.insert(info);
	}
	infile.close();
}

/**
 * readCyclesFile
 * @brief Read the cycles from a file.
 */
void NonAdjacentCycleStatsApp::readCyclesFile()
{
	std::ifstream infile;
	mCycles.clear();
	infile.open (mstrCyclesFile.c_str(), std::ifstream::in);
	std::string strLine;
	std::vector<std::string> fields;

	while(std::getline(infile, strLine).good())
	{
		annotate::cleanString(strLine, ' ');

		fields = annotate::splitStringFields(strLine, ":");

		std::string strPDBFile = fields[0];
		unsigned int uiModel = atol(fields[1].c_str());
		std::string strProfile = fields[2];
		std::string strPredProfile = fields[3];
		std::string strResIds = fields[4];
		std::string strSeq = fields[5];

		annotate::CycleProfile prof(strProfile);
		annotate::CycleProfile predProf(strPredProfile);

		annotate::CycleInfo::residue_profile resProfile;
		resProfile = getStrandResidues(strResIds, predProf);

		std::vector<std::string> residues = annotate::splitStringFields(strSeq, "-");
		annotate::CycleInfo cInfo(strPDBFile, uiModel, predProf, prof, resProfile, residues);

		mCycles.insert(cInfo);
	}

	infile.close();
}

std::vector<std::vector<std::string> > NonAdjacentCycleStatsApp::getStrandResidues(
	const std::string& aResidues,
	const annotate::CycleProfile& aProfile) const
{
	std::list<std::string> residues = getResidues(aResidues);
	std::vector<std::vector<std::string> > strandResidues;

	std::list<std::string>::const_iterator itRes = residues.begin();
	std::list<unsigned int>::const_iterator itProf;
	for(itProf = aProfile.strandProfile().begin();
		itProf != aProfile.strandProfile().end();
		++ itProf)
	{
		std::vector<std::string> strand;
		unsigned int iRes = 0;
		for(iRes = 0; iRes < *itProf; ++ iRes)
		{
			strand.push_back(*itRes);
			++ itRes;
		}
		strandResidues.push_back(strand);
		strand.clear();
	}
	return strandResidues;
}

std::list<std::string> NonAdjacentCycleStatsApp::getResidues(
	const std::string& aResidues) const
{
	std::list<std::string> residues;
	std::string strResidues = aResidues;

	annotate::cleanString(strResidues, ' ');
	while(0 < strResidues.size())
	{
		std::string strRes;
		std::size_t sep = strResidues.rfind('-');
		if(sep != std::string::npos)
		{
			strRes = strResidues.substr(sep + 1, strResidues.size() - (sep + 1));
			strResidues.erase(sep);
		}
		else
		{
			strRes = strResidues;
			strResidues.clear();
		}
		residues.push_front(strRes);
	}
	return residues;
}

std::set<annotate::CycleInfo> NonAdjacentCycleStatsApp::getNonAdjacentCycleFromModel(
	const GetModelRangeFunctor<annotate::InteractionInfo>::const_range& interactionRange,
	const GetModelRangeFunctor<annotate::CycleInfo>::const_range& cycleRange) const
{
	std::set<annotate::CycleInfo> returnVal;
	std::set<annotate::InteractionInfo>::const_iterator interIt;
	std::set<annotate::CycleInfo>::const_iterator cycleIt;
	for(interIt = interactionRange.first;
		interIt != interactionRange.second;
		++ interIt)
	{
		for(cycleIt = cycleRange.first; cycleIt != cycleRange.second; ++cycleIt)
		{
			if(cycleIt->contains(*interIt))
			{
				returnVal.insert(*cycleIt);
			}
		}
	}
	return returnVal;
}

std::set<annotate::CycleInfo> NonAdjacentCycleStatsApp::getNonAdjacentCycles() const
{
	GetModelRangeFunctor<annotate::CycleInfo> cycleRangeFunctor;
	GetModelRangeFunctor<annotate::InteractionInfo> interactionRangeFunctor;

	std::set<annotate::CycleInfo> returnVal;
	std::list<annotate::ModelInfo> models = this->getModels();

	std::list<annotate::ModelInfo>::const_iterator itModel;
	for(itModel = models.begin(); itModel != models.end(); ++itModel)
	{
		GetModelRangeFunctor<annotate::InteractionInfo>::const_range interRange;
		GetModelRangeFunctor<annotate::CycleInfo>::const_range cycleRange;

		interRange = interactionRangeFunctor(mInteractions, *itModel);
		cycleRange = cycleRangeFunctor(mCycles, *itModel);

		std::set<annotate::CycleInfo> cycles;
		cycles = getNonAdjacentCycleFromModel(interRange, cycleRange);

		returnVal.insert(cycles.begin(), cycles.end());
	}
	return returnVal;
}

// Get a list of models
std::list<annotate::ModelInfo> NonAdjacentCycleStatsApp::getModels() const
{
	const annotate::ModelInfo* pInfo = NULL;
	std::list<annotate::ModelInfo> models;
	std::set<annotate::InteractionInfo>::const_iterator it;
	for(it = mInteractions.begin(); it != mInteractions.end(); ++ it)
	{
		if(it != mInteractions.begin() && (*pInfo != it->getModelInfo()))
		{
			models.push_back(*pInfo);
		}
		pInfo = &(it->getModelInfo());
	}
	if(NULL != pInfo)
	{
		models.push_back(*pInfo);
	}
	return models;
}

std::string NonAdjacentCycleStatsApp::toString() const
{
	std::ostringstream oss;
	oss << "Number of interactions found : " << interactions().size() << std::endl;
	oss << "Number of cycle found : " << cycles().size() << std::endl;
	oss << "Number of non-adjacent cycle found : " << mNACycles.size() << std::endl;
	oss << "Number of cycles excluding non-adjacent ones : " << mSecondaryStructureCycles.size() << std::endl;
	oss << "Number of ignored cycle due to open strand " << muiOpenConnections << std::endl;
	oss << "Number of ignored cycle due to partial strand coverage " << muiPartialConnections << std::endl;
	oss << "Number of ignored cycle due to subcycles " << muiEnclosingCycles << std::endl;
	oss << "Number of ignored cycle due to multibranch connection " << muiMultibranchCycles << std::endl;
	oss << "Number of ignored cycle due to loose connection " << muiLooseCycles << std::endl;
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		oss << "------------------------------" << std::endl;
		oss << it->toString() << std::endl;
	}
	oss << statisticsToString() << std::endl;
	oss << cycleStatsToString() << std::endl;
	oss << mInteractionTable.toString() << std::endl;
	oss << interactionsToString() << std::endl;

	return oss.str();
}

std::set<annotate::CycleInfo> NonAdjacentCycleStatsApp::getCyclesWithInteraction(
	const GetModelRangeFunctor<annotate::CycleInfo>::const_range aRange,
	const std::set<annotate::Interaction>& aInteractions) const
{
	std::set<annotate::CycleInfo> cycles;
	std::set<annotate::CycleInfo>::const_iterator itCycle;
	for(itCycle = aRange.first; itCycle != aRange.second; ++ itCycle)
	{
		if(itCycle->shareInteraction(aInteractions))
		{
			cycles.insert(*itCycle);
		}
	}
	return cycles;
}

std::set<NACycleInfo> NonAdjacentCycleStatsApp::getConnectedCycles(
	const std::set<annotate::CycleInfo>& aNACycles) const
{
	GetModelRangeFunctor<annotate::CycleInfo> cycleRangeFunctor;
	std::set<NACycleInfo> oNACycles;
	std::set<annotate::CycleInfo>::const_iterator itCycle = aNACycles.begin();
	for(; itCycle != aNACycles.end(); ++ itCycle)
	{
		NACycleInfo cycle = *itCycle;
		unsigned int uiNbStrands = cycle.getNbStrands();
		for(unsigned int iStrand = 0; iStrand < uiNbStrands; ++ iStrand)
		{
			std::set<annotate::Interaction> interactions;
			interactions = cycle.getStrandInteractions(iStrand);
			GetModelRangeFunctor<annotate::CycleInfo>::const_range range;
			range = cycleRangeFunctor(mSecondaryStructureCycles, itCycle->getModelInfo());
			std::set<annotate::CycleInfo> connectedCycles;
			connectedCycles =  getCyclesWithInteraction(range, interactions);
			cycle.getStrandConnections(iStrand) = connectedCycles;
		}
		oNACycles.insert(cycle);
	}
	return oNACycles;
}

// Get interacting cycles pairs
std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > >
NonAdjacentCycleStatsApp::getDoubleInteractionCyclesPairsFromInteractionRange(
	const GetModelRangeFunctor<NAInteractionInfo>::const_range& aInterRange) const
{
	std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> > returnVal;
	std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> >  interactingCycles;
	std::set<NAInteractionInfo>::const_iterator it;
	for(it = aInterRange.first; it != aInterRange.second; ++ it)
	{
		// Get the interacting cycles from a single interaction
		std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> > modelCycles;
		modelCycles = getInteractingCyclesPairsFromInteraction(*it);

		// For each pair in the set, if it was already present (from another
		// interaction), add this pairing to the interacting cycles
		std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> >::const_iterator it2;
		for(it2 = modelCycles.begin(); it2 != modelCycles.end(); ++ it2)
		{
			std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> >::const_iterator it3;
			it3 = interactingCycles.find(*it2);
			if(it3 != interactingCycles.end())
			{
				returnVal.insert(*it3);
			}
			interactingCycles.insert(*it2);
		}
	}
	return returnVal;
}

// Get interacting cycles pairs for a model
std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > >
NonAdjacentCycleStatsApp::getDoubleInteractionCyclesPairsFromModel(
	const std::set<NAInteractionInfo>& aInteractions,
	const annotate::ModelInfo& aModel) const
{
	// Get the range of interactions from the model
	GetModelRangeFunctor<NAInteractionInfo> interactionRangeFunctor;
	GetModelRangeFunctor<NAInteractionInfo>::const_range interRange;
	std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> >  interactingCycles;
	interRange = interactionRangeFunctor(aInteractions, aModel);

	// Return the interactions for that range
	return getDoubleInteractionCyclesPairsFromInteractionRange(interRange);
}

// Get the set of all interacting cycles pairs with at least two interactions
// between them.
std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > > NonAdjacentCycleStatsApp::getDoubleInteractionCyclesPairs(
	const std::set<annotate::InteractionInfo>& aInteractions,
	const std::set<annotate::CycleInfo>& aCycles) const
{
	std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> > returnVal;
	std::set<NAInteractionInfo> interactions = getNAInteractions(aInteractions, aCycles);
	std::list<annotate::ModelInfo> models = this->getModels();
	std::list<annotate::ModelInfo>::const_iterator itModel;
	for(itModel = models.begin(); itModel != models.end(); ++itModel)
	{
		// Get the pair of cycles with at least two interactions in this model
		std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> > modelCycles;
		modelCycles = getDoubleInteractionCyclesPairsFromModel(
			interactions,
			*itModel);

		// Extract the cycles in the return value
		std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> >::const_iterator itPair;
		for(itPair = modelCycles.begin(); itPair != modelCycles.end(); ++ itPair)
		{
			returnVal.insert(*itPair);
		}
	}
	return returnVal;
}

// Get the set of all interacting cycles with at least two interactions between
// them
std::set<annotate::CycleInfo> NonAdjacentCycleStatsApp::getDoubleInteractionCycles() const
{
	std::set<annotate::CycleInfo> returnCycles;

	// Get the double interaction pairs
	std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > > pairs;
	pairs = getDoubleInteractionCyclesPairs(mInteractions, mSecondaryStructureCycles);

	// Extract the cycles in the return value
	std::set< std::pair<annotate::CycleInfo, annotate::CycleInfo> >::const_iterator itPair;
	for(itPair = pairs.begin(); itPair != pairs.end(); ++ itPair)
	{
		returnCycles.insert(itPair->first);
		returnCycles.insert(itPair->second);
	}

	return returnCycles;
}

std::set<CycleStatsEntry, lessCycleStatsEntry> NonAdjacentCycleStatsApp::computeCycleStatsInstances() const
{
	std::set<CycleStatsEntry, lessCycleStatsEntry> cycleCount;
	std::set<annotate::CycleInfo>::const_iterator itCycle;
	for(itCycle = mSecondaryStructureCycles.begin(); itCycle != mSecondaryStructureCycles.end(); ++ itCycle)
	{
		std::string strProfile = itCycle->getProfile().toString();
		std::string strSequence = itCycle->residuesString("");
		CycleStatsEntry statsEntry(strProfile, strSequence);
		std::set<CycleStatsEntry>::iterator itStats;
		itStats = cycleCount.find(statsEntry);
		if(itStats != cycleCount.end())
		{
			CycleStatsEntry entry = *itStats;
			cycleCount.erase(itStats);
			entry.addInstance(*itCycle);
			cycleCount.insert(entry);

		}else
		{
			CycleStatsEntry entry(strProfile, strSequence);
			entry.addInstance(*itCycle);
			cycleCount.insert(entry);
		}
	}
	return cycleCount;
}

void NonAdjacentCycleStatsApp::addInteractionBetween(std::set<CycleStatsEntry, lessCycleStatsEntry>& aCycleStats,
	const annotate::CycleInfo& aCycle1,
	const annotate::CycleInfo& aCycle2)
{
	std::string strProfile = aCycle1.getProfile().toString();
	std::string strSequence = aCycle1.residuesString("");
	CycleStatsEntry statsEntry(strProfile, strSequence);
	std::set<CycleStatsEntry>::iterator itStats;
	itStats = aCycleStats.find(statsEntry);
	if(itStats != aCycleStats.end())
	{
		CycleStatsEntry entry = *itStats;
		aCycleStats.erase(itStats);
		entry.addInteraction(aCycle1, aCycle2);
		mInteractionTable.addInteraction(aCycle1, aCycle2);
		aCycleStats.insert(entry);
	}
}

void NonAdjacentCycleStatsApp::computeCycleStatsInteractions(
	std::set<CycleStatsEntry, lessCycleStatsEntry>& aCycleStats,
	const std::set<annotate::InteractionInfo>& aInteractions,
	const std::set<annotate::CycleInfo>& aCycles)
{
	std::set<std::pair<annotate::CycleInfo, annotate::CycleInfo> > pairs;
	pairs = getDoubleInteractionCyclesPairs(aInteractions, aCycles);

	std::set<std::pair<annotate::CycleInfo, annotate::CycleInfo> >::const_iterator itPair;
	for(itPair = pairs.begin(); itPair != pairs.end(); ++ itPair)
	{
		addInteractionBetween(aCycleStats, itPair->first, itPair->second);
		addInteractionBetween(aCycleStats, itPair->second, itPair->first);
	}
}

// Get the interacting cycles from a single interaction
std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > >
NonAdjacentCycleStatsApp::getInteractingCyclesPairsFromInteraction(
	const NAInteractionInfo& aInteraction) const
{
	std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > > returnVal;
	std::set<annotate::CycleInfo>::const_iterator it1;
	std::set<annotate::CycleInfo>::const_iterator it2;
	it1 = aInteraction.fivePrimeCycles().begin();
	for(; it1 != aInteraction.fivePrimeCycles().end(); ++ it1)
	{
		it2 = aInteraction.threePrimeCycles().begin();
		for(; it2 != aInteraction.threePrimeCycles().end(); ++ it2)
		{
			std::pair< annotate::CycleInfo, annotate::CycleInfo > entry(*it1, *it2);
			if(*it2 < *it1)
			{
				entry.first = *it2;
				entry.second = *it1;
			}
			returnVal.insert(entry);
		}
	}
	return returnVal;
}

std::set<NAInteractionInfo> NonAdjacentCycleStatsApp::getNAInteractions(
	const std::set<annotate::InteractionInfo>& aInteractions,
	const std::set<annotate::CycleInfo>& aCycles) const
{
	GetModelRangeFunctor<annotate::CycleInfo> cycleRangeFunctor;
	GetModelRangeFunctor<annotate::InteractionInfo> interactionRangeFunctor;
	std::set<NAInteractionInfo> returnVal;
	std::list<annotate::ModelInfo> models = this->getModels();

	std::list<annotate::ModelInfo>::const_iterator itModel;
	for(itModel = models.begin(); itModel != models.end(); ++itModel)
	{
		GetModelRangeFunctor<annotate::InteractionInfo>::const_range interRange;
		GetModelRangeFunctor<annotate::CycleInfo>::const_range cycleRange;

		interRange = interactionRangeFunctor(aInteractions, *itModel);
		cycleRange = cycleRangeFunctor(aCycles, *itModel);

		std::set<NAInteractionInfo> cyclesAndPairs;
		cyclesAndPairs = getNAInteractionsFromModel(interRange, cycleRange);

		returnVal.insert(cyclesAndPairs.begin(), cyclesAndPairs.end());
	}
	return returnVal;
}

std::set<NAInteractionInfo> NonAdjacentCycleStatsApp::getNAInteractionsFromModel(
	const GetModelRangeFunctor<annotate::InteractionInfo>::const_range& interactionRange,
	const GetModelRangeFunctor<annotate::CycleInfo>::const_range& cycleRange) const
{
	std::set<NAInteractionInfo> returnVal;
	std::set<annotate::InteractionInfo>::const_iterator interIt;
	std::set<annotate::CycleInfo>::const_iterator cycleIt;
	for(interIt = interactionRange.first;
		interIt != interactionRange.second;
		++ interIt)
	{
		NAInteractionInfo oInfo(*interIt);
		for(cycleIt = cycleRange.first; cycleIt != cycleRange.second; ++cycleIt)
		{
			if(cycleIt->contains(interIt->getRes1()))
			{
				oInfo.addFivePrimeCycle(*cycleIt);
			} else if(cycleIt->contains(interIt->getRes2()))
			{
				oInfo.addThreePrimeCycle(*cycleIt);
			}
		}
		if((0 < oInfo.fivePrimeCycles().size()) && (0 < oInfo.threePrimeCycles().size()))
		{
			returnVal.insert(oInfo);
		}
	}
	return returnVal;
}

void NonAdjacentCycleStatsApp::computeCyclesStats()
{
	std::set<annotate::CycleInfo> interactingSecondaryStructureCycles;
	mCyclesStats = computeCycleStatsInstances();
	interactingSecondaryStructureCycles = this->getDoubleInteractionCycles();
	computeCycleStatsInteractions(mCyclesStats, mInteractions, interactingSecondaryStructureCycles);
}

void NonAdjacentCycleStatsApp::filterOutOpenConnections()
{
	unsigned int uiNbConnectedCycles = mConnectedCycles.size();
	bool bPass;
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		bPass = true;
		std::vector<std::set<annotate::CycleInfo> >::const_iterator itStrand;
		for(itStrand = it->getConnections().begin();
			itStrand != it->getConnections().end() && bPass;
			++itStrand)
		{
			if(0 == itStrand->size())
			{
				bPass = false;
			}
		}
		if(bPass)
		{
			filtered.insert(*it);
		}
	}
	mConnectedCycles = filtered;
	muiOpenConnections = uiNbConnectedCycles - mConnectedCycles.size();
}

void NonAdjacentCycleStatsApp::filterOutPartialCoverage()
{
	unsigned int uiNbConnectedCycles = mConnectedCycles.size();
	bool bPass;
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		bPass = true;
		NACycleInfo cycle(
			it->getPDBFile(),
			it->getModel(),
			it->getFileProfile(),
			it->getProfile(),
			it->getStrandResIds(),
			it->getSequence());

		unsigned int uiIndex;
		bool bPassStrand = true;
		for(uiIndex = 0;
			uiIndex < it->getConnections().size() && bPassStrand;
			++ uiIndex)
		{
			bPassStrand = false;
			std::set<annotate::Interaction> strand = it->getStrandInteractions(uiIndex);
			std::set<annotate::CycleInfo> connection = it->getConnections()[uiIndex];
			for(std::set<annotate::CycleInfo>::const_iterator itCycle = connection.begin();
				itCycle != connection.end();
				++ itCycle)
			{
				if(itCycle->hasStrandCoveringInteractions(strand))
				{
					cycle.getStrandConnections(uiIndex).insert(*itCycle);
					bPassStrand = true;
				}
			}
			bPass  = bPassStrand;
		}
		if(bPass)
		{
			filtered.insert(cycle);
		}
	}
	mConnectedCycles = filtered;
	muiPartialConnections = uiNbConnectedCycles - mConnectedCycles.size();
}

void NonAdjacentCycleStatsApp::filterOutEnclosingCycles()
{
	unsigned int uiNbConnectedCycles = mConnectedCycles.size();
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		NACycleInfo cycle(
			it->getPDBFile(),
			it->getModel(),
			it->getFileProfile(),
			it->getProfile(),
			it->getStrandResIds(),
			it->getSequence());

		unsigned int uiIndex;
		for(uiIndex = 0; uiIndex < it->getConnections().size(); ++ uiIndex)
		{
			std::set<annotate::CycleInfo> connection = it->getConnections()[uiIndex];
			for(std::set<annotate::CycleInfo>::const_iterator itCycle = connection.begin();
				itCycle != connection.end();
				++ itCycle)
			{
				bool bHasSubCycle = false;
				std::set<annotate::CycleInfo>::const_iterator itSubCycle;
				for(itSubCycle = connection.begin();
					itSubCycle != connection.end() && !bHasSubCycle;
					++ itSubCycle)
				{
					if(itCycle != itSubCycle)
					{
						if(itSubCycle->isSubCycleOf(*itCycle))
						{
							bHasSubCycle = itCycle->isSubCycleOf(*itSubCycle);
						}
					}
				}
				if(!bHasSubCycle)
				{
					cycle.getStrandConnections(uiIndex).insert(*itCycle);
				}
			}
		}
		filtered.insert(cycle);
	}
	mConnectedCycles = filtered;
	muiEnclosingCycles = uiNbConnectedCycles - mConnectedCycles.size();
}

void NonAdjacentCycleStatsApp::filterOutMultibranchCycles()
{
	unsigned int uiNbConnectedCycles = mConnectedCycles.size();
	bool bPass;
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		bPass = true;

		unsigned int uiIndex;
		for(uiIndex = 0;
			uiIndex < it->getConnections().size() && bPass;
			++ uiIndex)
		{
			std::set<annotate::Interaction> strand = it->getStrandInteractions(uiIndex);
			std::set<annotate::CycleInfo> connection = it->getConnections()[uiIndex];
			for(std::set<annotate::CycleInfo>::const_iterator itCycle = connection.begin();
				itCycle != connection.end() && bPass;
				++ itCycle)
			{
				if(itCycle->getProfile().type() == annotate::Cycle::eMULTIBRANCH)
				{
					bPass = false;
				}
			}
		}
		if(bPass)
		{
			filtered.insert(*it);
		}
	}
	mConnectedCycles = filtered;
	muiMultibranchCycles = uiNbConnectedCycles - mConnectedCycles.size();
}

void NonAdjacentCycleStatsApp::filterOutLooseCycles()
{
	unsigned int uiNbConnectedCycles = mConnectedCycles.size();
	bool bPass;
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		bPass = true;

		unsigned int uiIndex;
		for(uiIndex = 0;
			uiIndex < it->getConnections().size() && bPass;
			++ uiIndex)
		{
			std::set<annotate::Interaction> strand = it->getStrandInteractions(uiIndex);
			std::set<annotate::CycleInfo> connection = it->getConnections()[uiIndex];
			for(std::set<annotate::CycleInfo>::const_iterator itCycle = connection.begin();
				itCycle != connection.end() && bPass;
				++ itCycle)
			{
				if(itCycle->getProfile().type() == annotate::Cycle::eLOOSE)
				{
					bPass = false;
				}
			}
		}
		if(bPass)
		{
			filtered.insert(*it);
		}
	}
	mConnectedCycles = filtered;
	muiLooseCycles = uiNbConnectedCycles - mConnectedCycles.size();
}

void NonAdjacentCycleStatsApp::splitAdjacency()
{
	std::set<NACycleInfo> cycles;
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		assert(it->getConnections().size() == 1 || it->getConnections().size() == 2);
		if(1 == it->getConnections().size())
		{
			std::set<annotate::CycleInfo>::const_iterator itCycle;
			for(itCycle = it->getConnections()[0].begin();
				itCycle != it->getConnections()[0].end();
				++ itCycle)
			{
				NACycleInfo cycle(
					it->getPDBFile(),
					it->getModel(),
					it->getFileProfile(),
					it->getProfile(),
					it->getStrandResIds(),
					it->getSequence());
				cycle.getStrandConnections(0).insert(*itCycle);
				cycles.insert(cycle);
			}
		}
		else if(2 == it->getConnections().size())
		{
			std::set<annotate::CycleInfo>::const_iterator itCycle1;
			for(itCycle1 = it->getConnections()[0].begin();
				itCycle1 != it->getConnections()[0].end();
				++ itCycle1)
			{
				std::set<annotate::CycleInfo>::const_iterator itCycle2;
				for(itCycle2 = it->getConnections()[1].begin();
					itCycle2 != it->getConnections()[1].end();
					++ itCycle2)
				{
					NACycleInfo cycle(
						it->getPDBFile(),
						it->getModel(),
						it->getFileProfile(),
						it->getProfile(),
						it->getStrandResIds(),
						it->getSequence());
					cycle.getStrandConnections(0).insert(*itCycle1);
					cycle.getStrandConnections(1).insert(*itCycle2);
					cycles.insert(cycle);
				}
			}
		}
	}
	mConnectedCycles = cycles;
}

void NonAdjacentCycleStatsApp::compileStatistics()
{
	std::set<NACycleInfo>::const_iterator it;
	for(it = mConnectedCycles.begin(); it != mConnectedCycles.end(); ++ it)
	{
		std::string strEquation = it->equationString();
		std::map<std::string, unsigned int>::iterator it;
		it = mStatistics.find(strEquation);
		if(it == mStatistics.end())
		{
			mStatistics.insert(std::pair<std::string, unsigned int>(strEquation, 1));
		}
		else
		{
			it->second = it->second + 1;
		}
	}
}

std::string NonAdjacentCycleStatsApp::statisticsToString() const
{
	std::ostringstream oss;
	oss << "------------------------------" << std::endl;
	std::multimap<unsigned int, std::string> displayStats;
	std::map<std::string, unsigned int>::const_iterator itInput;
	for(itInput = mStatistics.begin(); itInput != mStatistics.end(); ++ itInput)
	{
		displayStats.insert(std::make_pair(itInput->second, itInput->first));
	}

	std::multimap<unsigned int, std::string>::const_iterator it;
	for(it = displayStats.begin(); it != displayStats.end(); ++ it)
	{
		oss << it->first << " : " << it->second << std::endl;
	}
	return oss.str();
}

// -----------------------------------------------------------------------------
std::string NonAdjacentCycleStatsApp::cycleStatsToString() const
{
	std::ostringstream oss;

	std::set<CycleStatsEntry>::const_iterator it;
	for(it = mCyclesStats.begin(); it != mCyclesStats.end(); ++it)
	{
		oss << it->toString() << std::endl;
	}
	return oss.str();
}

// -----------------------------------------------------------------------------
std::string NonAdjacentCycleStatsApp::interactionsToString() const
{
	std::ostringstream oss;

	std::list<std::string> profileList = mInteractionTable.getProfileList();
	std::list<std::string>::const_iterator it1 = profileList.begin();
	for(; it1 != profileList.end(); ++ it1)
	{
		std::list<std::string>::const_iterator it2;
		for(it2 = it1; it2 != profileList.end(); ++ it2)
		{
			oss << "//--------------------------------------------------------------" << std::endl;
			oss << *it1 << " vs " << *it2 << std::endl;
			oss << "//--------------------------------------------------------------" << std::endl;
			const InteractionTable::interacting_set& interacting = mInteractionTable(*it1, *it2);
			oss << interactionStats(interacting);
			InteractionTable::interacting_set::const_iterator it3;
			for(it3 = interacting.begin(); it3 != interacting.end(); ++ it3)
			{
				oss << it3->first.getModelInfo().getPDBFile() << " " << it3->first.getModelInfo().getModel() << " : ";
				oss << "(" << it3->first.toString() << ") vs (";
				oss << it3->second.toString() << ")" << std::endl;
			}
			if(!mstrOutputDirectory.empty())
			{
				outputInteractingCyclesFiles(*it1, *it2, interacting);
			}
		}
	}
	return oss.str();
}

std::string NonAdjacentCycleStatsApp::interactionStats(
	const InteractionTable::interacting_set& aInteracting) const
{
	std::ostringstream oss;
	std::map<std::pair<std::string, std::string>, unsigned int> interactions;
	InteractionTable::interacting_set::const_iterator it;
	for(it = aInteracting.begin(); it != aInteracting.end(); ++ it)
	{
		std::map<std::pair<std::string, std::string>, unsigned int>::iterator itMap;
		std::pair<std::string, std::string> interactingPair(it->first.residuesString(), it->second.residuesString());
		itMap = interactions.find(interactingPair);
		if(itMap == interactions.end())
		{
			interactions.insert(std::pair<std::pair<std::string, std::string>, unsigned int>(interactingPair, 1));
		}else
		{
			itMap->second += 1;
		}
	}
	std::map<std::pair<std::string, std::string>, unsigned int>::iterator itMap2;
	for(itMap2 = interactions.begin(); itMap2 != interactions.end(); ++ itMap2)
	{
		oss << itMap2->first.first << " => " << itMap2->first.second << " : ";
		oss << itMap2->second << std::endl;
	}
	return oss.str();
}

void NonAdjacentCycleStatsApp::outputInteractingCyclesFiles(
	const std::string& astrProfile1,
	const std::string& astrProfile2,
	const InteractionTable::interacting_set& aInteracting) const
{
	std::ostringstream oss;
	std::string strPrefix;
	InteractionTable::interacting_set::const_iterator it;
	std::ofstream outfile;
	oss << mstrOutputDirectory << "/" << astrProfile1 << "_vs_" << astrProfile2;
	oss << ".cycles";
	outfile.open (oss.str().c_str(), ios_base::out);
	for(it = aInteracting.begin(); it != aInteracting.end(); ++ it)
	{
		outfile << it->first.getPDBFile();
		outfile << " : " << it->first.getModelInfo().getModel();
		outfile << " : " << it->first.getFileProfile().toString();
		outfile << " : " << it->first.getProfile().toString();
		outfile << " : " << it->first.resIdsString();
		outfile << " : " << it->first.residuesString();
		outfile << std::endl;

		outfile << it->second.getPDBFile();
		outfile << " : " << it->second.getModelInfo().getModel();
		outfile << " : " << it->second.getFileProfile().toString();
		outfile << " : " << it->second.getProfile().toString();
		outfile << " : " << it->second.resIdsString();
		outfile << " : " << it->second.residuesString();
		outfile << std::endl;
	}
	outfile.close();
}

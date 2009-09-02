//                              -*- Mode: C++ -*-
// mcnacs.cc
// Copyright © 2009 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Jul 29 10:28:00 2009

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mcnacs.h"

#include "CycleInfo.h"
#include "Interaction.h"
#include "InteractionInfo.h"
#include "ModelInfo.h"
#include "NACycleInfo.h"

// libmcannotate
#include "AlgorithmExtra.h"
#include "CycleProfile.h"
#include "StringUtil.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>

bool gbSplitInteractions = false;
bool gbRemoveComposites = false;
std::string gstrCyclesFile;
std::string gstrSecondaryCyclesFile;
std::string gstrPairsFile;
const char* shortopts = "Vhus:p:c:";

typedef std::multimap<ModelInfo, InteractionInfo>::const_iterator inter_map_iterator;
typedef std::set<CycleInfo>::const_iterator cycle_set_iterator;
typedef std::pair<cycle_set_iterator, cycle_set_iterator> cycle_set_range;

// PROTOTYPES ------------------------------------------------------------------
void displayConnectedCycle(const NACycleInfo& aCycle);
std::string cycleInfoString(const CycleInfo& aCycle);

void version ()
{
	std::cout
		<< PACKAGE << " " << VERSION << " (" << __DATE__ << ")"
		<< std::endl;
}


void usage ()
{
	std::cout
	 	<< "usage: " << PACKAGE
		<< " [-huvV] -p <distant pairs file> -c <cycles file> [-s <secondary structure cycles file>]"
		<< std::endl;
}

void help ()
{
	std::cout
		<< "This program read cycle structures and return the corresponding residue ids."
		<< std::endl
		<< "  -h	print this help" << std::endl
		<< "  -u	split the non-adjacent cycles into sets of unique interacting cycle per strand" << std::endl
		<< "  -p	file containing the non-adjacent interacting pairs" << std::endl
		<< "  -c	file containing the identified cycles" << std::endl
		<< "  -s	file containing a second set of cycles identified in the secondary structure" << std::endl
		<< "  -V	print the software version info" << std::endl;
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
		case 'c':
		{
			gstrCyclesFile = optarg;
			break;
		}
		case 'f':
		{
			gbRemoveComposites = true;
			break;
		}
		case 'p':
		{
			gstrPairsFile = optarg;
			break;
		}
		case 's':
		{
			gstrSecondaryCyclesFile = optarg;
			break;
		}
		case 'u':
		{
			gbSplitInteractions = true;
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
}

std::string cycleResiduesString(const CycleInfo& aCycle)
{
	std::ostringstream oss;
	std::vector<std::string> resIds = aCycle.getResIds();
	std::vector<std::string>::const_iterator it = resIds.begin();
	for(it = resIds.begin(); it != resIds.end(); ++ it)
	{
		if(it != resIds.begin())
		{
			oss << "-";
		}
		oss << *it;
	}
	return oss.str();
}

std::multimap<ModelInfo, InteractionInfo> readPairsFile(const std::string& aFile)
{
	std::ifstream infile;
	std::multimap<ModelInfo, InteractionInfo> infos;

	infile.open (aFile.c_str(), std::ifstream::in);
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

		InteractionInfo info(strPDBFile, uiModel, strRes1, strRes2);
		infos.insert(std::pair<ModelInfo, InteractionInfo>(info.getModelInfo(), info));
	}

	infile.close();

	return infos;
}

std::list<std::string> getResidues(const std::string& aResidues)
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

std::vector<std::vector<std::string> > getStrandResidues(
	const std::string& aResidues,
	const annotate::CycleProfile& aProfile)
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

std::set<CycleInfo> readCyclesFile(const std::string& aFile)
{
	unsigned int i = 1;
	std::ifstream infile;
	std::set<CycleInfo> infos;
	infile.open (aFile.c_str(), std::ifstream::in);
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

		CycleInfo::residue_profile resProfile;
		resProfile = getStrandResidues(strResIds, prof);

		std::vector<std::string> residues = annotate::splitStringFields(strSeq, "-");
		CycleInfo cInfo(strPDBFile, uiModel, prof, resProfile, residues);

		infos.insert(cInfo);
		assert(infos.size() == i);
		++ i;
	}

	infile.close();

	return infos;
}

std::string modelInfoString(const CycleInfo& aCycle)
{
	std::ostringstream oss;
	oss << aCycle.getModelInfo().getPDBFile();
	oss << "[" << aCycle.getModelInfo().getModel() << "]";
	return oss.str();
}

std::string cycleInfoResiduesString(const CycleInfo& aCycle)
{
	std::ostringstream oss;

	std::vector<std::string> sequence = aCycle.getSequence();
	std::vector<std::string>::const_iterator it;
	for(it = sequence.begin(); it != sequence.end(); ++ it)
	{
		if(it != sequence.begin())
		{
			oss << "-";
		}
		oss << *it;
	}
	return oss.str();
}

std::string cycleInfoResIdsString(const CycleInfo& aCycle)
{
	std::ostringstream oss;

	CycleInfo::residue_profile::const_iterator itResStrand;
	for(itResStrand = aCycle.getStrandResIds().begin();
		itResStrand != aCycle.getStrandResIds().end();
		++ itResStrand)
	{
		if(itResStrand != aCycle.getStrandResIds().begin())
		{
			oss << ",";
		}
		oss << "{";
		CycleInfo::residue_strand::const_iterator itRes;
		for(itRes = itResStrand->begin();
			itRes != itResStrand->end();
			++ itRes)
		{
			if(itRes != itResStrand->begin())
			{
				oss << "-";
			}
			oss << *itRes;
		}
		oss << "}";
	}

	return oss.str();
}

std::string cycleInfoString(const CycleInfo& aCycle)
{
	std::ostringstream oss;

	annotate::CycleProfile profile = aCycle.getProfile();
	oss << profile.toString();
	oss << "\t: ";
	oss << cycleInfoResIdsString(aCycle);
	oss << " : ";
	oss << cycleInfoResiduesString(aCycle);
	return oss.str();
}

void displayCycleInfos(const std::set<CycleInfo>& aCycleInfos)
{
	std::set<CycleInfo>::const_iterator it = aCycleInfos.begin();
	for(it = aCycleInfos.begin(); it != aCycleInfos.end(); ++ it)
	{
		std::cout << modelInfoString(*it) << "\t" << cycleInfoString(*it);
		std::cout << std::endl;
	}
}

std::list<ModelInfo> getModels(
	const std::multimap<ModelInfo, InteractionInfo> aInfos)
{
	const ModelInfo* pInfo = NULL;
	std::list<ModelInfo> models;
	std::multimap<ModelInfo, InteractionInfo>::const_iterator it;
	for(it = aInfos.begin(); it != aInfos.end(); ++ it)
	{
		if(it != aInfos.begin() && (*pInfo != it->first))
		{
			models.push_back(*pInfo);
		}
		pInfo = &(it->first);
	}
	if(NULL != pInfo)
	{
		models.push_back(*pInfo);
	}
	return models;
}

std::set<CycleInfo> getNonAdjacentCycleFromModel(
	const std::pair<inter_map_iterator, inter_map_iterator>& interactionRange,
	const cycle_set_range& cycleRange)
{
	std::set<CycleInfo> oNACycles;
	inter_map_iterator interIt;
	cycle_set_iterator cycleIt;
	for(interIt = interactionRange.first;
		interIt != interactionRange.second;
		++ interIt)
	{
		for(cycleIt = cycleRange.first; cycleIt != cycleRange.second; ++cycleIt)
		{
			if(cycleIt->contains(interIt->second))
			{
				oNACycles.insert(*cycleIt);
			}
		}
	}
	return oNACycles;
}

std::pair<cycle_set_iterator, cycle_set_iterator> getModelRange(
	const std::set<CycleInfo>& aCycles,
	const ModelInfo& aModelInfo)
{
	std::pair<cycle_set_iterator, cycle_set_iterator> range(aCycles.end(), aCycles.end());
	cycle_set_iterator it = aCycles.begin();
	while(it != aCycles.end() && it->getModelInfo() < aModelInfo)
	{
		++ it;
	}
	range.first = it;
	while(it != aCycles.end() && !(aModelInfo < it->getModelInfo()))
	{
		++ it;
	}
	range.second = it;
	return range;
}

std::set<CycleInfo> getNonAdjacentCycle(
	const std::multimap<ModelInfo, InteractionInfo>& aInteractions,
	const std::set<CycleInfo>& aCycles)
{
	std::set<CycleInfo> oNACycles;

	std::list<ModelInfo> models = getModels(aInteractions);

	std::list<ModelInfo>::const_iterator itModel;
	for(itModel = models.begin(); itModel != models.end(); ++itModel)
	{
		std::pair<inter_map_iterator, inter_map_iterator> interRange;
		std::pair<cycle_set_iterator, cycle_set_iterator> cycleRange;

		interRange = aInteractions.equal_range(*itModel);
		cycleRange = getModelRange(aCycles, *itModel);

		std::set<CycleInfo> oModelCycle;
		oModelCycle = getNonAdjacentCycleFromModel(interRange, cycleRange);

		oNACycles.insert(oModelCycle.begin(), oModelCycle.end());
	}
	return oNACycles;
}

std::set<CycleInfo> removeCycles(
	std::multimap<ModelInfo, CycleInfo>& aCycles,
	const std::multimap<ModelInfo, CycleInfo>& aToRemove)
{
	std::set<CycleInfo> cycles;
	std::multimap<ModelInfo, CycleInfo>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		cycles.insert(it->second);
	}

	std::set<CycleInfo> cyclesToRemove;
	for(it = aToRemove.begin(); it != aToRemove.end(); ++ it)
	{
		cyclesToRemove.insert(it->second);
	}

	return annotate::SetDifference(cycles, cyclesToRemove);
}

std::set<CycleInfo> getCyclesWithInteraction(
	const cycle_set_range aRange,
	const std::set<Interaction>& aInteractions)
{
	std::set<CycleInfo> cycles;
	std::set<CycleInfo>::const_iterator itCycle;
	for(itCycle = aRange.first; itCycle != aRange.second; ++ itCycle)
	{
		if(itCycle->shareInteraction(aInteractions))
		{
			cycles.insert(*itCycle);
		}
	}
	return cycles;
}

std::set<NACycleInfo> getConnectedCycles(
	const std::set<CycleInfo>& aNACycles,
	const std::set<CycleInfo>& aACycles)
{
	std::set<NACycleInfo> oNACycles;
	std::set<CycleInfo>::const_iterator itCycle = aNACycles.begin();
	for(; itCycle != aNACycles.end(); ++ itCycle)
	{
		NACycleInfo cycle = *itCycle;
		unsigned int uiNbStrands = cycle.getNbStrands();
		for(unsigned int iStrand = 0; iStrand < uiNbStrands; ++ iStrand)
		{
			std::set<Interaction> interactions;
			interactions = cycle.getStrandInteractions(iStrand);
			std::pair<cycle_set_iterator, cycle_set_iterator> range;
			range = getModelRange(aACycles, itCycle->getModelInfo());
			std::set<CycleInfo> connectedCycles;
			connectedCycles =  getCyclesWithInteraction(range, interactions);
			cycle.getStrandConnections(iStrand) = connectedCycles;
		}
		oNACycles.insert(cycle);
	}
	return oNACycles;
}

std::string connectedCycleEquationString(const NACycleInfo& aCycle)
{
	std::ostringstream oss;
	std::string strProfile = aCycle.getProfile().toString();
	oss << strProfile << " = ";
	std::vector< std::set<CycleInfo> >::const_iterator it;
	for(it = aCycle.getConnections().begin();
		it != aCycle.getConnections().end();
		++ it)
	{
		if(it != aCycle.getConnections().begin())
		{
			oss << " + ";
		}
		if(1 < it->size())
		{
			oss << "(";
		}
		std::set<CycleInfo>::const_iterator itCycle;
		for(itCycle = it->begin(); itCycle != it->end(); ++ itCycle)
		{
			if(itCycle != it->begin())
			{
				oss << " + ";
			}
			oss << itCycle->getProfile().toString();
		}
		if(1 < it->size())
		{
			oss << ")";
		}
	}
	return oss.str();
}


std::set<NACycleInfo> filterOutOpenConnections(std::set<NACycleInfo>& aCycles)
{
	bool bPass;
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		bPass = true;
		std::vector< std::set<CycleInfo> >::const_iterator itStrand;
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
	return filtered;
}

std::set<NACycleInfo> filterOutPartialCoverage(
	const std::set<NACycleInfo>& aCycles)
{
	bool bPass;
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		bPass = true;
		NACycleInfo cycle(
			it->getPDBFile(),
			it->getModel(),
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
			std::set<Interaction> strand = it->getStrandInteractions(uiIndex);
			std::set<CycleInfo> connection = it->getConnections()[uiIndex];
			for(std::set<CycleInfo>::const_iterator itCycle = connection.begin();
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
	return filtered;
}

std::set<NACycleInfo> filterOutEnclosingCycles(
	const std::set<NACycleInfo>& aCycles)
{
	std::set<NACycleInfo> filtered;
	std::set<NACycleInfo>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		NACycleInfo cycle(
			it->getPDBFile(),
			it->getModel(),
			it->getProfile(),
			it->getStrandResIds(),
			it->getSequence());

		unsigned int uiIndex;
		for(uiIndex = 0; uiIndex < it->getConnections().size(); ++ uiIndex)
		{
			std::set<CycleInfo> connection = it->getConnections()[uiIndex];
			for(std::set<CycleInfo>::const_iterator itCycle = connection.begin();
				itCycle != connection.end();
				++ itCycle)
			{
				bool bHasSubCycle = false;
				std::set<CycleInfo>::const_iterator itSubCycle;
				for(itSubCycle = connection.begin();
					itSubCycle != connection.end() && !bHasSubCycle;
					++ itSubCycle)
				{
					if(itCycle != itSubCycle)
					{
						bHasSubCycle = itSubCycle->isSubCycleOf(*itCycle);
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
	return filtered;
}

void displayConnectedCycle(const NACycleInfo& aCycle)
{
	std::cout << "index:" << std::endl;
	std::cout << "\t" << connectedCycleEquationString(aCycle) << std::endl;
	std::cout << "info:" << std::endl;
	std::cout << "\t" << modelInfoString(aCycle) << std::endl;
	std::cout << "\t" << cycleInfoString(aCycle) << std::endl;
	std::cout << "\t=" << std::endl;
	std::vector< std::set<CycleInfo> >::const_iterator it;
	for(it = aCycle.getConnections().begin();
		it != aCycle.getConnections().end();
		++ it)
	{
		if(it != aCycle.getConnections().begin())
		{
			std::cout << "\t+" << std::endl;
		}
		std::set<CycleInfo>::const_iterator itCycle;
		for(itCycle = it->begin(); itCycle != it->end(); ++itCycle)
		{
			std::cout << "\t" << cycleInfoString(*itCycle) << std::endl;
		}
	}
}

void displayConnectedCycles(const std::set<NACycleInfo>& aCycles)
{
	std::set<NACycleInfo>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		std::cout << "------------------------------" << std::endl;
		displayConnectedCycle(*it);
	}
}

std::set<NACycleInfo> splitAdjacency(const std::set<NACycleInfo>& aCycles)
{
	std::set<NACycleInfo> cycles;
	std::set<NACycleInfo>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		assert(it->getConnections().size() == 1 || it->getConnections().size() == 2);
		if(1 == it->getConnections().size())
		{
			std::set<CycleInfo>::const_iterator itCycle;
			for(itCycle = it->getConnections()[0].begin();
				itCycle != it->getConnections()[0].end();
				++ itCycle)
			{
				NACycleInfo cycle(
					it->getPDBFile(),
					it->getModel(),
					it->getProfile(),
					it->getStrandResIds(),
					it->getSequence());
				cycle.getStrandConnections(0).insert(*itCycle);
				cycles.insert(cycle);
			}
		}
		else if(2 == it->getConnections().size())
		{
			std::set<CycleInfo>::const_iterator itCycle1;
			for(itCycle1 = it->getConnections()[0].begin();
				itCycle1 != it->getConnections()[0].end();
				++ itCycle1)
			{
				std::set<CycleInfo>::const_iterator itCycle2;
				for(itCycle2 = it->getConnections()[1].begin();
					itCycle2 != it->getConnections()[1].end();
					++ itCycle2)
				{
					NACycleInfo cycle(
						it->getPDBFile(),
						it->getModel(),
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
	return cycles;
}

std::map<std::string, unsigned int> compileStatistics(
	const std::set<NACycleInfo>& aCycles)
{
	std::map<std::string, unsigned int> stats;
	std::set<NACycleInfo>::const_iterator it;
	for(it = aCycles.begin(); it != aCycles.end(); ++ it)
	{
		std::string strEquation = connectedCycleEquationString(*it);
		std::map<std::string, unsigned int>::iterator it;
		it = stats.find(strEquation);
		if(it == stats.end())
		{
			stats.insert(std::pair<std::string, unsigned int>(strEquation, 1));
		}
		else
		{
			it->second = it->second + 1;
		}
	}
	return stats;
}

void displayStatistics(const std::map<std::string, unsigned int>& aStats)
{
	std::cout << "------------------------------" << std::endl;
	std::multimap<unsigned int, std::string> displayStats;
	std::map<std::string, unsigned int>::const_iterator itInput;
	for(itInput = aStats.begin(); itInput != aStats.end(); ++ itInput)
	{
		displayStats.insert(std::make_pair(itInput->second, itInput->first));
	}

	std::multimap<unsigned int, std::string>::const_iterator it;
	for(it = displayStats.begin(); it != displayStats.end(); ++ it)
	{
		std::cout << it->first << " : " << it->second << std::endl;
	}
}

int main (int argc, char *argv[])
{
	read_options (argc, argv);

	std::multimap<ModelInfo, InteractionInfo> interactionInfos = readPairsFile(gstrPairsFile);
	std::set<CycleInfo> cycleInfos = readCyclesFile(gstrCyclesFile);
	std::set<CycleInfo> cycleNAInfos = getNonAdjacentCycle(interactionInfos, cycleInfos);

	std::cout << "Number of interactions found : " << interactionInfos.size() << std::endl;
	std::cout << "Number of cycle found : " << cycleInfos.size() << std::endl;
	if(!gstrSecondaryCyclesFile.empty())
	{
		std::set<CycleInfo> secondInfo = readCyclesFile(gstrSecondaryCyclesFile);
		std::cout << "Number of supplementary cycle found : " << secondInfo.size() << std::endl;

		cycleInfos = annotate::SetUnion(cycleInfos, secondInfo);
		std::cout << "Total number of potential cycle found : " << cycleInfos.size() << std::endl;
	}
	std::cout << "Number of non-adjacent cycle found : " << cycleNAInfos.size() << std::endl;

	// Remove the non-adjacent cycle
	cycleInfos = annotate::SetDifference(cycleInfos, cycleNAInfos);
	std::cout << "Number of cycles excluding non-adjacent ones : " << cycleInfos.size() << std::endl;

	std::set<NACycleInfo> connectedCycles;
	connectedCycles = getConnectedCycles(cycleNAInfos, cycleInfos);

	// Remove the cycles not interacting with other NCMs
	unsigned int uiNbConnectedCycles = connectedCycles.size();
	connectedCycles = filterOutOpenConnections(connectedCycles);
	std::cout << "Number of ignored cycle due to open strand ";
	std::cout << (uiNbConnectedCycles - connectedCycles.size()) << std::endl;

	// Remove the cycles not covering entirely the relation
	uiNbConnectedCycles = connectedCycles.size();
	connectedCycles = filterOutPartialCoverage(connectedCycles);
	std::cout << "Number of ignored cycle due to partial strand coverage ";
	std::cout << (uiNbConnectedCycles - connectedCycles.size()) << std::endl;

	// Remove the cycles who has subcycles in the same interactions
	// ( e.g. NGNRAN vs GNRA, keep GNRA only )
	uiNbConnectedCycles = connectedCycles.size();
	connectedCycles = filterOutEnclosingCycles(connectedCycles);
	std::cout << "Number of ignored cycle due to subcycles ";
	std::cout << (uiNbConnectedCycles - connectedCycles.size()) << std::endl;

	if(gbSplitInteractions)
	{
		connectedCycles = splitAdjacency(connectedCycles);
	}

	displayConnectedCycles(connectedCycles);

	std::map<std::string, unsigned int> stats = compileStatistics(connectedCycles);
	displayStatistics(stats);

	return EXIT_SUCCESS;
}

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
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>

#include "CycleInfo.h"
#include "InteractionInfo.h"
#include "ModelInfo.h"

std::string gstrCyclesFile;
std::string gstrPairsFile;
const char* shortopts = "Vhp:c:";

typedef std::multimap<ModelInfo, InteractionInfo>::const_iterator inter_map_iterator;
typedef std::multimap<ModelInfo, CycleInfo>::const_iterator cycle_map_iterator;


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
		<< " [-hlvV] -p <distant pairs file> -c <cycles file>"
		<< std::endl;
}

void help ()
{
	std::cout	
		<< "This program read cycle structures and return the corresponding residue ids." 
		<< std::endl
		<< "  -h	print this help" << std::endl
		<< "  -p	file containing the non-adjacent interacting pairs" << std::endl
		<< "  -c	file containing the identified cycles" << std::endl
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
		case 'p':
		{
			gstrPairsFile = optarg;
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

std::multimap<ModelInfo, InteractionInfo> readPairsFile(const std::string& aFile)
{
	std::ifstream infile;
	std::multimap<ModelInfo, InteractionInfo> infos;

	infile.open (aFile.c_str(), std::ifstream::in);
	std::string strLine;

	while (std::getline(infile, strLine).good())
	{
		cleanString(strLine, ' ');
		
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

std::list<unsigned int> getProfile(const std::string& astrProfile)
{
	std::list<unsigned int> profile;
	std::string::const_iterator it;
	for(it = astrProfile.begin(); it != astrProfile.end(); ++ it)
	{
		std::string strNumber(1, *it);
		profile.push_back(atol (strNumber.c_str()));
	}
	return profile;	
}

std::list<std::string> getResidues(const std::string& aResidues)
{
	std::list<std::string> residues;
	std::string strResidues = aResidues;
	
	cleanString(strResidues, ' ');
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
	const std::list<unsigned int>& aProfile)
{
	std::list<std::string> residues = getResidues(aResidues);
	std::vector<std::vector<std::string> > strandResidues;
	
	std::list<std::string>::const_iterator itRes = residues.begin();
	std::list<unsigned int>::const_iterator itProf;
	for(itProf = aProfile.begin(); itProf != aProfile.end(); ++ itProf)
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

std::multimap<ModelInfo, CycleInfo> readCyclesFile(const std::string& aFile)
{
	std::ifstream infile;
	std::multimap<ModelInfo, CycleInfo> infos;
	infile.open (aFile.c_str(), std::ifstream::in);
	std::string strLine;
	
	while(std::getline(infile, strLine).good())
	{
		cleanString(strLine, ' ');
		
		std::size_t sep = strLine.rfind(':');
		std::string strSeq = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);
		
		sep = strLine.rfind(':');
		std::string strResIds = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);
		
		sep = strLine.rfind(':');
		std::string strProfile = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);
		
		sep = strLine.rfind(':');
		std::string strPredProfile = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);
		
		sep = strLine.rfind(':');
		std::string strModel = strLine.substr(sep + 1, strLine.size() - (sep + 1));
		strLine.erase(sep);
		unsigned int uiModel = atol(strModel.c_str());
		
		std::string strPDBFile = strLine;
		std::list<unsigned int> prof = getProfile(strProfile);
		CycleInfo::residue_profile resProfile;
		resProfile = getStrandResidues(strResIds, prof); 
		CycleInfo cInfo(strPDBFile, uiModel, resProfile);
		infos.insert(std::pair<ModelInfo, CycleInfo>(cInfo.getModelInfo(), cInfo));
	}
	
	infile.close();
	
	return infos;
}

void displayCycleInfos(const std::multimap<ModelInfo, CycleInfo>& aCycleInfos)
{
	std::multimap<ModelInfo, CycleInfo>::const_iterator it = aCycleInfos.begin();
	for(it = aCycleInfos.begin(); it != aCycleInfos.end(); ++ it)
	{
		const CycleInfo* pInfo = &(it->second);
		std::cout	<< pInfo->getPDBFile() << ":"
					<< pInfo->getModel() << ":";
		
		std::vector<unsigned int> profile = pInfo->getProfile();
		std::vector<unsigned int>::const_iterator itProf;
		for(itProf = profile.begin(); itProf != profile.end(); ++ itProf)
		{
			std::cout << *itProf;
		}
		std::cout << ":";
		CycleInfo::residue_profile::const_iterator itResStrand;
		for(itResStrand = pInfo->getStrandResidues().begin(); 
			itResStrand != pInfo->getStrandResidues().end(); 
			++ itResStrand)
		{
			if(itResStrand != pInfo->getStrandResidues().begin())
			{
				std::cout << ",";
			}
			std::cout << "{";
			CycleInfo::residue_strand::const_iterator itRes;
			for(itRes = itResStrand->begin(); 
				itRes != itResStrand->end(); 
				++ itRes)
			{
				if(itRes != itResStrand->begin())
				{
					std::cout << "-";
				}
				std::cout << *itRes;
			}
			std::cout << "}";
		}
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
	return models;
}

std::multimap<ModelInfo, CycleInfo> getNonAdjacentCycleFromModel(
	const std::pair<inter_map_iterator, inter_map_iterator>& interactionRange,
	const std::pair<cycle_map_iterator, cycle_map_iterator>& cycleRange)
{
	std::multimap<ModelInfo, CycleInfo> oNACycles;
	inter_map_iterator interIt;
	cycle_map_iterator cycleIt;
	for(interIt = interactionRange.first; 
		interIt != interactionRange.second; 
		++ interIt) 
	{
		for(cycleIt = cycleRange.first; cycleIt != cycleRange.second; ++cycleIt)
		{
			if(cycleIt->second.contains(interIt->second))
			{
				oNACycles.insert(*cycleIt);				
			}
		}
	}
	return oNACycles;
}

std::multimap<ModelInfo, CycleInfo> getNonAdjacentCycle(
	const std::multimap<ModelInfo, InteractionInfo>& aInteractions,
	const std::multimap<ModelInfo, CycleInfo>& aCycles)
{
	std::multimap<ModelInfo, CycleInfo> oNACycles;
	
	std::list<ModelInfo> models = getModels(aInteractions);
	
	std::list<ModelInfo>::const_iterator itModel;
	for(itModel = models.begin(); itModel != models.end(); ++itModel)
	{
		std::pair<inter_map_iterator, inter_map_iterator> interRange;
		std::pair<cycle_map_iterator, cycle_map_iterator> cycleRange;
		
		interRange = aInteractions.equal_range(*itModel);
		cycleRange = aCycles.equal_range(*itModel);
		
		std::multimap<ModelInfo, CycleInfo> oModelCycle;
		oModelCycle = getNonAdjacentCycleFromModel(interRange, cycleRange);
		
		oNACycles.insert(oModelCycle.begin(), oModelCycle.end());		
	}
	return oNACycles;
}
	
int main (int argc, char *argv[])
{
	read_options (argc, argv);
	
	std::multimap<ModelInfo, InteractionInfo> Interactioninfos = readPairsFile(gstrPairsFile);
	std::multimap<ModelInfo, CycleInfo> cycleInfos = readCyclesFile(gstrCyclesFile);
	
	std::multimap<ModelInfo, CycleInfo> cycleNAInfos = getNonAdjacentCycle(
		Interactioninfos, 
		cycleInfos);
	
	displayCycleInfos(cycleNAInfos);
	
	return EXIT_SUCCESS;	
}

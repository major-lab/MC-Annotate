/*
 * RNAGroupFile.cc
 *
 *  Created on: Apr 22, 2010
 *      Author: blanchmf
 */

#include "RNAGroupFile.h"

#include "StringUtil.h"

#include "CycleInfo.h"

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>

namespace annotate {

RNAGroupFile::RNAGroupFile()
{}

RNAGroupFile::~RNAGroupFile()
{
	mGroups.clear();
}

void RNAGroupFile::read(const std::string& astrFileName)
{
	std::ifstream infile;
	infile.open(astrFileName.c_str(), std::ios_base::in);
	mGroups.clear();
	if(infile.good())
	{
		std::string strLine;
		while(std::getline(infile, strLine).good())
		{
			readLine(strLine);
		}

		// Renumber offsets
		renumberOffsets();
	}
	else
	{
		// TODO : This should be an exception
		std::cout << "Error opening file " << astrFileName << std::endl;
	}
	infile.close();
}

/**
 * readLine
 * @brief Read a line from an RNAGroup file.
  */
void RNAGroupFile::readLine(const std::string& astrLine)
{
	std::string strLine = astrLine;
	cleanString(strLine, ' ');
	std::vector<std::string> fields = splitStringFields(strLine, ":");
	unsigned int uiGroup = atol(fields[0].c_str());
	int iOffset = atoi(fields[1].c_str());

	// Read the files
	for(unsigned int i = 2; i < fields.size(); ++ i)
	{
		std::vector<std::string> fileFields = splitStringFields(fields[i], "_");
		std::string strPDB = fileFields[0];
		unsigned int uiModel = atol(fileFields[1].c_str());
		char cChain = fileFields[2][0];
		RNAGroupFileEntry entry(
			strPDB,
			uiModel,
			cChain,
			iOffset);
		mGroups[uiGroup].push_back(entry);
		mEntryGroupMap[entry] = uiGroup;
	}
}

void RNAGroupFile::renumberOffsets()
{
	std::map<unsigned int, std::list<RNAGroupFileEntry> >::iterator it;
	for(it = mGroups.begin(); it != mGroups.end(); ++ it)
	{
		int iMinOffset = 0;
		std::list<RNAGroupFileEntry>::iterator itFile;

		// Find the minimum offset
		for(itFile = it->second.begin(); itFile != it->second.end(); ++ itFile)
		{
			if(itFile == it->second.begin())
			{
				iMinOffset = itFile->offset();
			}
			else
			{
				iMinOffset = std::min(itFile->offset(), iMinOffset);
			}
		}

		// Renumber so all offsets are >= 0
		for(itFile = it->second.begin(); itFile != it->second.end(); ++ itFile)
		{
			itFile->offset() = itFile->offset() - iMinOffset;
		}
	}
}

unsigned int RNAGroupFile::getGroup(const RNAGroupFileEntry& aEntry) const
{
	std::map<RNAGroupFileEntry, unsigned int>::const_iterator it;
	it = mEntryGroupMap.find(aEntry);
	assert(it != mEntryGroupMap.end()); // TODO : Make this an exception
	return it->second;
}

unsigned int RNAGroupFile::getGroup(const CycleInfo& aCycle) const
{
	std::map<RNAGroupFileEntry, unsigned int>::const_iterator it;
	std::set<char> chains = aCycle.getChains();
	assert(1 == chains.size()); // TODO : Make this an exception
	RNAGroupFileEntry entry(aCycle.getPDBFile(), aCycle.getModel(), *(chains.begin()), 0);
	it = mEntryGroupMap.find(entry);
	assert(it != mEntryGroupMap.end()); // TODO : Make this an exception
	return it->second;
}

unsigned int RNAGroupFile::getGroup(const InteractionInfo& aInteraction) const
{
	std::map<RNAGroupFileEntry, unsigned int>::const_iterator it;
	std::set<char> chains = aInteraction.getChains();
	assert(1 == chains.size()); // TODO : Make this an exception
	RNAGroupFileEntry entry(aInteraction.getModelInfo().getPDBFile(), aInteraction.getModelInfo().getModel(), *(chains.begin()), 0);
	it = mEntryGroupMap.find(entry);
	assert(it != mEntryGroupMap.end()); // TODO : Make this an exception
	return it->second;
}

};


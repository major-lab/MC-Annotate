//                              -*- Mode: C++ -*-
// t4buildtable.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Sun May 2 12:13:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t4tablebuilderh1.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

T4BuildTableH1::T4BuildTableH1(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
: T4TableBuilder(astrGroupFile, astrInteractions)
{
}

void T4BuildTableH1::computeStatistics()
{
	// Initialize the count and frequencies to 0
	initializeTables();

	// Count the number of interactions of each type
	countInteractions();

	// Compute the frequencies
	computeFrequencies();
}

void T4BuildTableH1::initializeTables()
{
	// Initialize the count and frequencies to 0
	mInteractionsCount.resize(4);
	mFrequencies.resize(4);
	for(unsigned int i = 0; i < 4; ++ i)
	{
		mInteractionsCount[i].resize(4);
		mFrequencies[i].resize(4);
		for(unsigned int j = 0; j < 4; ++ j)
		{
			mInteractionsCount[i][j].resize(mInteractionTypeMap.size(), 0);
			mFrequencies[i][j].resize(mInteractionTypeMap.size(), 0.0f);
		}
	}
}

void T4BuildTableH1::countInteractions()
{
	// Count the number of interactions of each type
	std::set<interaction_cycle_pair>::const_iterator it;
	for(it = mInteractions.begin(); it != mInteractions.end(); ++ it)
	{
		unsigned int uiGroup = mGroupFile.getGroup(it->first);
		assert(uiGroup == mGroupFile.getGroup(it->second));
		assert(uiGroup == mGroupFile.getGroup(it->third));
		if(mGroupModelMap[uiGroup].isNew(*it))
		{
			unsigned int uiTypeId = interactionTypeId(it->first);
			unsigned int uiFlipTypeId = flipInteractionTypeId(it->first);
			unsigned int uiNucId1 = getNucleotideId(it->first.resId1(), it->second);
			unsigned int uiNucId2 = getNucleotideId(it->first.resId2(), it->third);
			mInteractionsCount[uiNucId1][uiNucId2][uiTypeId] += 1;
			mInteractionsCount[uiNucId2][uiNucId1][uiFlipTypeId] += 1;
		}
	}
}

void T4BuildTableH1::computeFrequencies()
{
	for(unsigned int i = 0; i < 4; ++ i)
	{
		for(unsigned int j = 0; j < 4; ++ j)
		{
			unsigned int uiCount = 0;
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				uiCount += mInteractionsCount[i][j][k];
			}
			if(0 < uiCount)
			{
				// Count the number of a given type for the given pair
				for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
				{
					mFrequencies[i][j][k] = (float)(mInteractionsCount[i][j][k]) / (float)(uiCount);
				}
			}
		}
	}
}

unsigned int T4BuildTableH1::getNucleotideId(
	const mccore::ResId& aResId,
	const annotate::CycleInfo& aCycle) const
{
	unsigned int uiNucId = 0;
	std::string strNuc = aCycle.getNucleotideString(aResId);
	if(strNuc == "A")
	{
		uiNucId = 0;
	}else if(strNuc == "C")
	{
		uiNucId = 1;
	}
	else if(strNuc == "G")
	{
		uiNucId = 2;
	}
	else if(strNuc == "U")
	{
		uiNucId = 3;
	}else
	{
		std::cout << "Invalid nucleotide : " << aResId << "," << strNuc << std::endl;
		std::cout << aCycle.residuesString("-");
		assert(false);
	}
	return uiNucId;
}

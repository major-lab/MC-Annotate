//                              -*- Mode: C++ -*-
// t1buildtable.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Wed Mar 17 11:14:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t2tablebuilderh2.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

T2TableBuilderH2::T2TableBuilderH2(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
: T2TableBuilder(astrGroupFile, astrInteractions)
{
}

void T2TableBuilderH2::computeStatistics()
{
	// Initialize the count and frequencies to 0
	initializeTables();

	// Count the number of interactions of each type
	countInteractions();

	// Compute the frequencies
	computeFrequencies();
}

void T2TableBuilderH2::initializeTables()
{
	// Initialize the count and frequencies to 0
	mInteractionsCount.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		mInteractionsCount[i].resize(mInteractionTypeMap.size(), 0);
	}

	mFrequencies.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		mInteractionsCount[i].resize(mProfileMap.size());
		mFrequencies[i].resize(mProfileMap.size());
		for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
		{

			mFrequencies[i][j].resize(mInteractionTypeMap.size(), 0.0f);
		}
	}
}

void T2TableBuilderH2::countInteractions()
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
			unsigned int uiProfile1 = profileId(it->second);
			unsigned int uiProfile2 = profileId(it->third);

			mInteractionsCount[uiProfile1][uiTypeId] += 1;
			mInteractionsCount[uiProfile2][uiFlipTypeId] += 1;
		}
	}
}

void T2TableBuilderH2::computeFrequencies()
{
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		// Count the number of a given type for the given pair
		unsigned int uiCountLeft = 0;

		for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
		{
			uiCountLeft += mInteractionsCount[i][k];
		}
		if(0 < uiCountLeft)
		{
			for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
			{
				unsigned int uiCountRight = 0;
				for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
				{
					uiCountRight += mInteractionsCount[j][k];
				}

				if(0 < uiCountRight)
				{
					for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
					{
						float fLeft = ((float)mInteractionsCount[i][k]) / ((float)uiCountLeft);
						float fRight = ((float)mInteractionsCount[j][k]) / ((float)uiCountRight);;
						mFrequencies[i][j][k] = std::min(fLeft, fRight);
					}
				}
			}
		}
	}
}


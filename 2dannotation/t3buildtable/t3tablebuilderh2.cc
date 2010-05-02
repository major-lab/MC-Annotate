//                              -*- Mode: C++ -*-
// t3buildtable.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Sat May 1 14:38:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t3tablebuilderh2.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

T3TableBuilderH2::T3TableBuilderH2(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
: T3TableBuilder(astrGroupFile, astrInteractions)
{
}

void T3TableBuilderH2::computeStatistics()
{
	// Initialize the count and frequencies to 0
	initializeTables();

	// Count the number of interactions of each type
	countInteractions();

	// Compute the frequencies
	computeFrequencies();
}

void T3TableBuilderH2::initializeTables()
{
	// Initialize the count to 0
	mInteractionsCount.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		unsigned int uiSize = profileSize(i);
		mInteractionsCount[i].resize(mInteractionTypeMap.size());
		for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
		{
			mInteractionsCount[i][k].resize(uiSize, 0);
		}
	}

	// Initialize the frequencies to 0
	mFrequencies.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		unsigned int uiSize1 = profileSize(i);
		mFrequencies[i].resize(mProfileMap.size());
		for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
		{
			unsigned int uiSize2 = profileSize(j);
			mFrequencies[i][j].resize(mInteractionTypeMap.size());
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				mFrequencies[i][j][k].resize(uiSize1);
				for(unsigned int l = 0; l < uiSize1; ++ l)
				{
					mFrequencies[i][j][k][l].resize(uiSize2, 0.0f);
				}
			}
		}
	}
}

void T3TableBuilderH2::countInteractions()
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
			interaction_coordinate coords = getInteractionCoordinates(*it);

			mInteractionsCount[uiProfile1][uiTypeId][coords.first] += 1;
			mInteractionsCount[uiProfile2][uiFlipTypeId][coords.second] += 1;
		}
	}
}

void T3TableBuilderH2::computeFrequencies()
{
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		// Count the number of a given type for the given pair
		unsigned int uiSize1 = profileSize(i);
		std::vector<unsigned int> countLeft;
		countLeft.resize(uiSize1, 0);
		for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
		{
			for(unsigned int l = 0; l < uiSize1; ++ l)
			{
				countLeft[l] += mInteractionsCount[i][k][l];
			}
		}
		for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
		{
			unsigned int uiSize2 = profileSize(j);
			std::vector<unsigned int> countRight;
			countRight.resize(uiSize2, 0);
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				for(unsigned int m = 0; m < uiSize2; ++ m)
				{
					countRight[m] += mInteractionsCount[j][k][m];
				}
			}

			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				for(unsigned int l = 0; l < uiSize1; ++ l)
				{
					if(0 < countLeft[l])
					{
						float fLeft = ((float)mInteractionsCount[i][k][l]) / ((float)countLeft[l]);
						for(unsigned int m = 0; m < uiSize2; ++ m)
						{
							if(0 < countRight[m])
							{
								float fRight = ((float)mInteractionsCount[j][k][m]) / ((float)countRight[m]);
								mFrequencies[i][j][k][l][m] = std::min(fLeft, fRight);
							}
						}
					}
				}
			}
		}
	}
}


//                              -*- Mode: C++ -*-
// t3buildtable.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Sat May 1 14:36:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t3tablebuilderh1.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

T3BuildTableH1::T3BuildTableH1(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
: T3TableBuilder(astrGroupFile, astrInteractions)
{
}

void T3BuildTableH1::computeStatistics()
{
	// Initialize the count and frequencies to 0
	initializeTables();

	// Count the number of interactions of each type
	countInteractions();

	// Compute the frequencies
	computeFrequencies();
}

void T3BuildTableH1::initializeTables()
{
	// Initialize the count and frequencies to 0
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

	// Initialize the count and frequencies to 0
	mInteractionsCount.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		mInteractionsCount[i].resize(mProfileMap.size());
		for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
		{
			unsigned int uiSize2 = profileSize(j);
			mInteractionsCount[i][j].resize(mInteractionTypeMap.size());
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				unsigned int uiSize1 = 0;
				std::pair<unsigned int, unsigned int> sizes = profileSizes(i);
				if(sizes.first == sizes.second)
				{
					uiSize1 = sizes.first;
				}else
				{
					uiSize1 = sizes.first + sizes.second;
				}
				mInteractionsCount[i][j][k].resize(uiSize1);
				for(unsigned int l = 0; l < mInteractionsCount[i][j][k].size(); ++ l)
				{
					std::pair<unsigned int, unsigned int> sizes2 = profileSizes(j);
					if(sizes2.first == sizes2.second)
					{
						uiSize2 = sizes2.first;
					}else
					{
						uiSize2 = sizes2.first + sizes2.second;
					}
					mInteractionsCount[i][j][k][l].resize(uiSize2, 0);
				}
			}
		}
	}
}

void T3BuildTableH1::countInteractions()
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
			mInteractionsCount[uiProfile1][uiProfile2][uiTypeId][coords.first][coords.second] += 1;
			mInteractionsCount[uiProfile2][uiProfile1][uiFlipTypeId][coords.second][coords.first] += 1;
		}
	}
}

void T3BuildTableH1::computeFrequencies()
{
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
		{
			// Count the number of a given type for the given pair
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				computeFrequencies(i, j, k);
			}
		}
	}
}

void T3BuildTableH1::computeFrequencies(
	unsigned int auiProfile1,
	unsigned int auiProfile2,
	unsigned int auiTypeId )
{
	unsigned int uiProfileSize1 = profileSize(auiProfile1);
	unsigned int uiProfileSize2 = profileSize(auiProfile2);

	// Count the number of interactions of each type at coordinates
	std::vector<std::vector<unsigned int> > count;
	count.resize(uiProfileSize1);
	for(unsigned int l = 0; l < uiProfileSize1; ++ l)
	{
		count[l].resize(uiProfileSize2, 0);
	}

	for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
	{
		for(unsigned int l = 0; l < uiProfileSize1; ++ l)
		{
			unsigned int uiLeftPos = getNucleotideSymetricPosition(auiProfile1, l);
			for(unsigned int m = 0; m < uiProfileSize2; ++ m)
			{
				unsigned int uiRightPos = getNucleotideSymetricPosition(auiProfile2, m);
				count[uiLeftPos][uiRightPos] += mInteractionsCount[auiProfile1][auiProfile2][k][uiLeftPos][uiRightPos]; // TODO : Make this unique
			}
		}
	}

	for(unsigned int l = 0; l < uiProfileSize1; ++ l)
	{
		unsigned int uiLeftPos = getNucleotideSymetricPosition(auiProfile1, l);
		for(unsigned int m = 0; m < uiProfileSize2; ++ m)
		{
			unsigned int uiRightPos = getNucleotideSymetricPosition(auiProfile2, m);
			if(0 < count[uiLeftPos][uiRightPos])
			{
				float fCountPos = (float)mInteractionsCount[auiProfile1][auiProfile2][auiTypeId][uiLeftPos][uiRightPos];
				float fFreq =  fCountPos / ((float)count[uiLeftPos][uiRightPos]);
				mFrequencies[auiProfile1][auiProfile2][auiTypeId][l][m] = fFreq;
			}
		}
	}
}



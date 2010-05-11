//                              -*- Mode: C++ -*-
// t3buildtableh3.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Sat May 1 12:20:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t3tablebuilderh3.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

T3TableBuilderH3::T3TableBuilderH3(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
: T3TableBuilder(astrGroupFile, astrInteractions)
{
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eWatson, annotate::BasePair::eCis)] = 0;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eWatson, annotate::BasePair::eTrans)] = 1;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eWatson, annotate::BasePair::eUnknown)] = 2;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eHoogsteen, annotate::BasePair::eCis)] = 3;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eHoogsteen, annotate::BasePair::eTrans)] = 4;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eHoogsteen, annotate::BasePair::eUnknown)] = 5;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eSugar, annotate::BasePair::eCis)] = 6;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eSugar, annotate::BasePair::eTrans)] = 7;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eSugar, annotate::BasePair::eUnknown)] = 8;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::eRibose, annotate::BasePair::eUnknown)] = 9;
	mOrientedFaceMap[oriented_face(annotate::InteractionInfo::ePhosphate, annotate::BasePair::eUnknown)] = 10;
}

void T3TableBuilderH3::computeStatistics()
{
	// Initialize the count and frequencies to 0
	initializeTables();

	// Count the number of interactions of each type
	countInteractions();

	// Compute the frequencies
	computeFrequencies();
}

void T3TableBuilderH3::initializeTables()
{
	// Initialize the count and frequencies to 0
	mInteractionsCount.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		std::pair<unsigned int, unsigned int> sizes = profileSizes(i);
		unsigned int uiSize = sizes.first;
		if(sizes.first != sizes.second)
		{
			uiSize += sizes.second;
		}
		mInteractionsCount[i].resize(mOrientedFaceMap.size());
		for(unsigned int k = 0; k < mOrientedFaceMap.size(); ++ k)
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

void T3TableBuilderH3::countInteractions()
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
			std::pair<unsigned int, unsigned int> faces;
			faces = getOrientedFacesIds(it->first);
			unsigned int uiProfile1 = profileId(it->second);
			unsigned int uiProfile2 = profileId(it->third);
			interaction_coordinate coords = getInteractionCoordinates(*it);

			mInteractionsCount[uiProfile1][faces.first][coords.first] += 1;
			mInteractionsCount[uiProfile2][faces.second][coords.second] += 1;
		}
	}
}

void T3TableBuilderH3::computeFrequencies()
{
	std::vector<std::vector<unsigned int> > count = computeInteractionsByPositions();
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		unsigned int uiProfile1 = profileSize(i);
		for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
		{
			unsigned int uiProfile2 = profileSize(j);
			for(unsigned int k = 0; k < mInteractionTypeMap.size(); ++ k)
			{
				for(unsigned int l = 0; l < uiProfile1; ++ l)
				{
					unsigned int uiLeftPos = getNucleotideSymetricPosition(i, l);
					unsigned int uiCount1 = count[i][uiLeftPos];
					if(0 < uiCount1)
					{
						for(unsigned int m = 0; m < uiProfile2; ++ m)
						{
							unsigned int uiRightPos = getNucleotideSymetricPosition(j, m);
							unsigned int uiCount2 = count[j][uiRightPos];
							if(0 < uiCount2)
							{
								std::pair<unsigned int, unsigned int> faces = getFacesIdFromTypeId(k);
								float fLeftCount = (float)mInteractionsCount[i][faces.first][uiLeftPos];
								float fRightCount = (float)mInteractionsCount[j][faces.second][uiRightPos];
								float fLeft = fLeftCount / ((float)uiCount1);
								float fRight = fRightCount / ((float)uiCount2);;
								mFrequencies[i][j][k][l][m] = std::min(fLeft, fRight);
							}
						}
					}
				}
			}
		}
	}
}

std::pair<unsigned int, unsigned int> T3TableBuilderH3::getOrientedFacesIds(
	const annotate::InteractionInfo& aInteraction) const
{
	oriented_face face1(aInteraction.face1(), aInteraction.orientation());
	oriented_face face2(aInteraction.face2(), aInteraction.orientation());
	return internalGetOrientedFacesIds(face1, face2);
}

std::pair<unsigned int, unsigned int> T3TableBuilderH3::internalGetOrientedFacesIds(
	const oriented_face& aFace1, const oriented_face& aFace2) const
{
	unsigned int uiTypeId1 = 0;
	unsigned int uiTypeId2 = 0;
	std::map<oriented_face, unsigned int>::const_iterator it;

	// Get the first face
	it = mOrientedFaceMap.find(aFace1);
	if(it != mOrientedFaceMap.end())
	{
		uiTypeId1 = it->second;
	}else
	{
		std::cout << "Could not find oriented face type : (" << aFace1.first;
				std::cout << "," << aFace1.second;
				std::cout << ")" << std::endl;
		assert(false);
	}

	// Get the second face
	it = mOrientedFaceMap.find(aFace2);
	if(it != mOrientedFaceMap.end())
	{
		uiTypeId2 = it->second;
	}else
	{
		std::cout << "Could not find oriented face type : (" << aFace2.first;
				std::cout << "," << aFace2.second;
				std::cout << ")" << std::endl;
		assert(false);
	}
	return std::pair<unsigned int, unsigned int>(uiTypeId1, uiTypeId2);
}

std::pair<unsigned int, unsigned int> T3TableBuilderH3::getFacesIdFromTypeId(unsigned int auiTypeId) const
{
	std::pair<unsigned int, unsigned int> facesId;
	bool bFound = false;
	std::map<interaction_type, unsigned int>::const_iterator it = mInteractionTypeMap.begin();
	while(it != mInteractionTypeMap.end() && !bFound )
	{
		if(auiTypeId == it->second)
		{
			interaction_type inter = it->first;
			face_pair facePair = inter.first;
			face eFace1 = facePair.first;
			face eFace2 = facePair.second;
			annotate::BasePair::enOrientation eOrientation = inter.second;
			oriented_face face1(eFace1, eOrientation);
			oriented_face face2(eFace2, eOrientation);
			facesId = internalGetOrientedFacesIds(face1, face2);
			bFound = true;
		}
		else
		{
			++ it;
		}
	}
	assert(it != mInteractionTypeMap.end());
	return facesId;
}

std::vector<std::vector<unsigned int> > T3TableBuilderH3::computeInteractionsByPositions() const
{
	std::vector<std::vector<unsigned int> > count;
	count.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		count[i].resize(mOrientedFaceMap.size(), 0);
		for(unsigned int k = 0; k < mOrientedFaceMap.size(); ++ k)
		{
			for(unsigned int l = 0; l < mInteractionsCount[i][k].size(); ++ l)
			{
				count[i][l] += mInteractionsCount[i][k][l];
			}
		}
	}
	return count;
}


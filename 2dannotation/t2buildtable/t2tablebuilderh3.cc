//                              -*- Mode: C++ -*-
// t1buildtableh3.cc
// Copyright © 2001-10 Laboratoire de Biologie Informatique et Théorique.
//                     Université de Montréal
// Author           : Marc-Frédérick Blanchet
// Created On       : Sat May 1 12:20:00 2010


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "t2tablebuilderh3.h"

#include "mccore/Messagestream.h"
#include "mccore/Version.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

T2TableBuilderH3::T2TableBuilderH3(
	const std::string& astrGroupFile,
	const std::string& astrInteractions)
: T2TableBuilder(astrGroupFile, astrInteractions)
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

void T2TableBuilderH3::computeStatistics()
{
	// Initialize the count and frequencies to 0
	initializeTables();

	// Count the number of interactions of each type
	countInteractions();

	// Compute the frequencies
	computeFrequencies();
}

void T2TableBuilderH3::initializeTables()
{
	// Initialize the count and frequencies to 0
	mInteractionsCount.resize(mProfileMap.size());
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		mInteractionsCount[i].resize(mOrientedFaceMap.size(), 0);
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

void T2TableBuilderH3::countInteractions()
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

			mInteractionsCount[uiProfile1][faces.first] += 1;
			mInteractionsCount[uiProfile2][faces.second] += 1;
		}
	}
}

void T2TableBuilderH3::computeFrequencies()
{
	for(unsigned int i = 0; i < mProfileMap.size(); ++ i)
	{
		// Count the number of a given type for the given pair
		unsigned int uiCountLeft = 0;

		for(unsigned int k = 0; k < mOrientedFaceMap.size(); ++ k)
		{
			uiCountLeft += mInteractionsCount[i][k];
		}
		if(0 < uiCountLeft)
		{
			for(unsigned int j = 0; j < mProfileMap.size(); ++ j)
			{
				unsigned int uiCountRight = 0;
				for(unsigned int k = 0; k < mOrientedFaceMap.size(); ++ k)
				{
					uiCountRight += mInteractionsCount[j][k];
				}

				if(0 < uiCountRight)
				{
					for(unsigned int l = 0; l < mInteractionTypeMap.size(); ++ l)
					{
						std::pair<unsigned int, unsigned int> faces = getFacesIdFromTypeId(l);
						float fLeft = ((float)mInteractionsCount[i][faces.first]) / ((float)uiCountLeft);
						float fRight = ((float)mInteractionsCount[j][faces.second]) / ((float)uiCountRight);;
						mFrequencies[i][j][l] = std::min(fLeft, fRight);
					}
				}
			}
		}
	}
}

std::pair<unsigned int, unsigned int> T2TableBuilderH3::getOrientedFacesIds(
	const annotate::InteractionInfo& aInteraction) const
{
	oriented_face face1(aInteraction.face1(), aInteraction.orientation());
	oriented_face face2(aInteraction.face2(), aInteraction.orientation());
	return internalGetOrientedFacesIds(face1, face2);
}

std::pair<unsigned int, unsigned int> T2TableBuilderH3::internalGetOrientedFacesIds(
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

std::pair<unsigned int, unsigned int> T2TableBuilderH3::getFacesIdFromTypeId(unsigned int auiTypeId) const
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


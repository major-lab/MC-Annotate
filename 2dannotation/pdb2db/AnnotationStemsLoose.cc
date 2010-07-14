/*
 * AnnotationStemsLoose.cc
 *
 *  Created on: Jul 5, 2010
 *      Author: blanchmf
 */

#include "AlgorithmExtra.h"
#include "AnnotationStemsLoose.h"
#include "AnnotateModel.h"
#include "AnnotationInteractions.h"

#include <cassert>
#include <sstream>

// Static members
std::string AnnotationStemsLoose::mstrAnnotationName = "StemsLoose";

	// Methods
AnnotationStemsLoose::AnnotationStemsLoose()
{
	addRequirement<annotate::AnnotationInteractions>();
}

AnnotationStemsLoose::~AnnotationStemsLoose()
{
	clear();
}

void AnnotationStemsLoose::update(annotate::AnnotateModel& aModel)
{
	// Remove previous annotation if any
	clear();

	// Create all the potential stems
	std::vector<annotate::Stem> stems;
	getPotentialStems(aModel, stems);

	// Keep the remaining stems for the annotation
	std::vector<annotate::Stem>::iterator itStem = stems.begin();
	for(itStem = stems.begin(); itStem != stems.end(); ++itStem)
	{
		// Name the stem properly
		std::string name = stemName(aModel, *itStem);
		itStem->name(name);
		mStems.push_back(*itStem);
	}
}

std::string AnnotationStemsLoose::stemName(
	const annotate::AnnotateModel& aModel,
	const annotate::Stem& aStem)
{
	std::ostringstream oss;
	mccore::ResId r1 = aStem.basePairs().front().fResId;
	mccore::ResId r2 = aStem.basePairs().back().fResId;
	mccore::ResId r3= aStem.basePairs().front().rResId;
	mccore::ResId r4 = aStem.basePairs().back().rResId;
	oss << r1 << "-" << r2 << ", " << r3 << "-" << r4;

	return oss.str();
}

std::string AnnotationStemsLoose::output() const
{
	std::ostringstream oss;
	vector< annotate::Stem >::const_iterator it;

	for (it = mStems.begin (); mStems.end () != it; ++it)
	{
		oss << it->name() << std::endl;
	}
	return oss.str();
}

const std::vector< annotate::Stem >& AnnotationStemsLoose::getStems() const
{
	return mStems;
}

void AnnotationStemsLoose::clear()
{
	mStems.clear();
}

void AnnotationStemsLoose::getPotentialStems(
	const annotate::AnnotateModel& aModel,
	std::vector<annotate::Stem>& aStems) const
{
	std::vector<annotate::Stem> potentialStems;

	// Get all the pairs
	std::vector< annotate::BasePair > allPairs;
	const annotate::AnnotationInteractions* pInteractions = NULL;
	pInteractions = aModel.getAnnotation<annotate::AnnotationInteractions>();
	assert(NULL != pInteractions);
	allPairs = pInteractions->pairs();

	// Filter out invalid pairs
	std::set<annotate::BasePair> pairs = filterPairs(allPairs);

	// Filter out multi-chain pairs
	pairs = filterOutMultiChainsPairs(pairs);

	std::map<mccore::ResId, std::set<annotate::BasePair> > potentialPairs;
	std::set<annotate::BasePair>::const_iterator itPair = pairs.begin();
	for(; itPair != pairs.end(); ++ itPair)
	{
		assert(itPair->fResId < itPair->rResId);
		potentialPairs[itPair->fResId].insert(*itPair);
	}

	std::map<mccore::ResId, std::set<annotate::BasePair> >::iterator itPotential;
	std::map<mccore::ResId, std::set<annotate::BasePair> >::iterator itNext;

	// Extract a first potential
	itPotential = potentialPairs.begin();
	while(itPotential != potentialPairs.end())
	{
		if(0 == itPotential->second.size())
		{
			potentialPairs.erase(itPotential);
			itPotential = potentialPairs.begin();
		}
		else
		{
			// Pick a base pair
			annotate::BasePair pair = *itPotential->second.begin();

			// Create a stem
			annotate::Stem stem;
			stem.push_back(pair);

			itNext = itPotential;
			itNext ++;
			bool bContinues = true;
			while(bContinues && itNext != potentialPairs.end())
			{
				bContinues = false;
				std::set<annotate::BasePair>::iterator itPair;
				for(itPair = itNext->second.begin(); itPair != itNext->second.end(); ++ itPair)
				{
					if(stem.continues(*itPair))
					{
						stem.push_back(*itPair);
						itNext->second.erase(itPair);
						itNext ++;
						bContinues = true;
						break;
					}
				}
			}

			// Keep the stem
			potentialStems.push_back(stem);

			// Remove the base pair
			itPotential->second.erase(pair);
			if(0 == itPotential->second.size())
			{
				potentialPairs.erase(itPotential);
			}

			// Redo the process
			itPotential = potentialPairs.begin();
		}
	}
	aStems = potentialStems;
}

bool AnnotationStemsLoose::shouldFilterPair(const annotate::BasePair& aPair) const
{
	bool bShouldFilter = true;
	annotate::BasePair::face_vector::const_iterator itFace;
	annotate::BasePair::face_vector faces = aPair.faces();
	for(itFace = faces.begin(); itFace != faces.end() && bShouldFilter; ++ itFace)
	{
		bool bFace1Ok = itFace->first->isW() || itFace->first->isH() || itFace->first->isS();
		bool bFace2Ok = itFace->second->isW() || itFace->second->isH() || itFace->second->isS();

		bShouldFilter = !(bFace1Ok && bFace2Ok);
	}
	return bShouldFilter;
}

std::set<annotate::BasePair> AnnotationStemsLoose::filterPairs(const std::vector<annotate::BasePair>& aPairs) const
{
	std::set< annotate::BasePair > filteredPairs;
	std::vector<annotate::BasePair>::const_iterator it = aPairs.begin();
	for(; it != aPairs.end(); ++ it)
	{
		bool bShouldFilter = shouldFilterPair(*it);
		if(!bShouldFilter)
		{
			filteredPairs.insert(*it);
		}
	}
	return filteredPairs;
}

std::set< annotate::BasePair > AnnotationStemsLoose::filterOutMultiChainsPairs(
	std::set< annotate::BasePair >& aPairs) const
{
	std::set< annotate::BasePair > filteredSet;
	for(std::set<annotate::BasePair>::const_iterator it = aPairs.begin();
		aPairs.end() != it;
		++ it)
	{
		if(it->fResId.getChainId() == it->rResId.getChainId())
		{
			filteredSet.insert(*it);
		}
	}
	return filteredSet;
}

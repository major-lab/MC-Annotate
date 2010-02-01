/*
 * NAInteractionInfo.cc
 *
 *  Created on: Nov 12, 2009
 *      Author: blanchmf
 */

#include "NAInteractionInfo.h"
#include <sstream>

NAInteractionInfo::NAInteractionInfo(const annotate::InteractionInfo& aInteraction)
: annotate::InteractionInfo(aInteraction)
{
}

void NAInteractionInfo::addFivePrimeCycle(const annotate::CycleInfo& aCycle)
{
	mFivePrimeCycles.insert(aCycle);
}

void NAInteractionInfo::addThreePrimeCycle(const annotate::CycleInfo& aCycle)
{
	mThreePrimeCycles.insert(aCycle);
}

std::string NAInteractionInfo::toString() const
{
	std::ostringstream ss;
	ss << "Interaction : (" << getRes1() << "," << getRes2() << ")" << std::endl;
	std::set<annotate::CycleInfo>::const_iterator it;
	ss << "5 prime " << mFivePrimeCycles.size() << std::endl;
	for(it = mFivePrimeCycles.begin(); it != mFivePrimeCycles.end(); ++it)
	{
		ss << it->toString() << std::endl;
	}
	ss << "3 prime " << mThreePrimeCycles.size() << std::endl;
	for(it = mThreePrimeCycles.begin(); it != mThreePrimeCycles.end(); ++it)
	{
		ss << it->toString() << std::endl;
	}
	return ss.str();
}

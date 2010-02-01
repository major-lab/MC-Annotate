/*
 * CycleStatsEntry.cc
 *
 *  Created on: Dec 16, 2009
 *      Author: blanchmf
 */

#include "CycleStatsEntry.h"
#include "CycleInfo.h"
#include <sstream>

CycleStatsEntry::CycleStatsEntry(
	const std::string& astrProfile,
	const std::string& astrSequence)
{
	mstrProfile = astrProfile;
	mstrSequence = astrSequence;
	muiInteractions = 0;
}

CycleStatsEntry::~CycleStatsEntry()
{
	mInstances.clear();
}

void CycleStatsEntry::addInstance(const annotate::CycleInfo& aCycle)
{
	std::pair<annotate::CycleInfo, cyclestatsentry_set> aPair;
	aPair.first = aCycle;
	mInstances.insert(aPair);
}

// Add an interaction between a cycle of this profile and another cycle
void CycleStatsEntry::addInteraction(
	const annotate::CycleInfo& aThisCycle,
	const annotate::CycleInfo& aOtherCycle)
{
	// Update the interaction count
	muiInteractions ++;

	// Get the instance entry
	std::map<annotate::CycleInfo, cyclestatsentry_set >::iterator it;
	it = mInstances.find(aThisCycle);

	if(it == mInstances.end())
	{
		std::pair<std::map<annotate::CycleInfo, cyclestatsentry_set>::iterator, bool> insertResult;
		std::pair<annotate::CycleInfo, cyclestatsentry_set> aPair;
		aPair.first = aThisCycle;
		insertResult = mInstances.insert(aPair);
		it = insertResult.first;
	}

	// Compute the stats on specific types of interactions
	std::string strProfile = aOtherCycle.getProfile().toString();
	std::string strSequence = aOtherCycle.residuesString("");
	CycleStatsEntry statsEntry(strProfile, strSequence);
	std::set<CycleStatsEntry>::iterator itStats;
	itStats = it->second.find(statsEntry);
	if(itStats != it->second.end())
	{
		CycleStatsEntry entry = *itStats;
		it->second.erase(itStats);
		entry.addInstance(aOtherCycle);
		it->second.insert(entry);
	}else
	{
		CycleStatsEntry entry(strProfile, strSequence);
		entry.addInstance(aOtherCycle);
		it->second.insert(entry);
	}
}

std::string CycleStatsEntry::toString() const
{
	std::ostringstream oss;
	oss << mstrProfile << " ";
	oss << mstrSequence << " ";
	oss << mInstances.size() << " ";
	oss << countInteractingInstances() << " ";
	oss << countInteractions() << " ";
	oss << "{";
	std::map<annotate::CycleInfo, cyclestatsentry_set >::const_iterator it;
	for(it = mInstances.begin(); it != mInstances.end(); ++it)
	{
		if(0 < it->second.size())
		{
			if(it != mInstances.begin()) oss << ", ";
			oss << "[(" << it->first.toString() << ") ";
			cyclestatsentry_set::const_iterator itEntry;
			for(itEntry = it->second.begin(); itEntry != it->second.end(); ++ itEntry)
			{
				if(itEntry != it->second.begin()) oss << ", ";
				oss << itEntry->mstrProfile << " " << itEntry->mstrSequence << " " << itEntry->mInstances.size();
			}
			oss << "]";
		}
	}
	oss << "}";
	return oss.str();
}

// Count the number of instances of this profile involved in interactions
unsigned int CycleStatsEntry::countInteractingInstances() const
{
	unsigned int uiCount = 0;
	std::map<annotate::CycleInfo, cyclestatsentry_set >::const_iterator it;
	for(it = mInstances.begin(); it != mInstances.end(); ++it)
	{
		if(0 < it->second.size())
		{
			uiCount ++;
		}
	}
	return uiCount;
}

// Count the number of interactions by instances of this profile
unsigned int CycleStatsEntry::countInteractions() const
{
	unsigned int uiCount = 0;
	std::map<annotate::CycleInfo, cyclestatsentry_set >::const_iterator it;
	for(it = mInstances.begin(); it != mInstances.end(); ++it)
	{
		uiCount += it->second.size();
	}
	return uiCount;
}

// Operator to compare two CycleStatsEntry
bool lessCycleStatsEntry::operator()(
	const CycleStatsEntry& aLeft,
	const CycleStatsEntry& aRight) const
{
	bool bLess = (aLeft.profile() < aRight.profile());
	bLess = bLess || ((aLeft.profile() == aRight.profile()) && (aLeft.sequence() < aRight.sequence()));
	return bLess;
}

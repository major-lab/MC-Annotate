/*
 * RNAGroupModel.cc
 *
 *  Created on: Apr 23, 2010
 *      Author: blanchmf
 */

#include "RNAGroupModel.h"

#include <cassert>

namespace annotate {

RNAGroupModel::RNAGroupModel(const std::list<RNAGroupFileEntry>& aFiles)
{
	std::list<RNAGroupFileEntry>::const_iterator it;
	for(it = aFiles.begin(); it != aFiles.end(); ++ it)
	{
		mFiles.insert(*it);
	}
}

bool RNAGroupModel::isNew(const CycleInfo& aCycle) const
{
	CycleInfo cycle = getGenericCycle(aCycle);
	std::set<CycleInfo>::const_iterator it = mCycles.find(cycle);
	return it == mCycles.end();
}

void RNAGroupModel::addCycle(const CycleInfo& aCycle)
{
	assert(internalIsNew(aCycle));
	CycleInfo cycle = getGenericCycle(aCycle);
	mCycles.insert(cycle);
}

RNAGroupFileEntry RNAGroupModel::getGroupFileEntry(
	const CycleInfo& aCycle) const
{
	std::set<char> chains = aCycle.getChains();
	assert(1 == chains.size());

	RNAGroupFileEntry entry(
		aCycle.getPDBFile(),
		aCycle.getModel(),
		*(chains.begin()),
		0);
	return entry;
}

RNAGroupFileEntry RNAGroupModel::getGroupFileEntry(
	const InteractionInfo& aInter) const
{
	std::set<char> chains = aInter.getChains();
	assert(1 == chains.size());

	RNAGroupFileEntry entry(
		aInter.getModelInfo().getPDBFile(),
		aInter.getModelInfo().getModel(),
		*(chains.begin()),
		0);
	return entry;
}

CycleInfo RNAGroupModel::getGenericCycle(const CycleInfo& aCycle) const
{
	RNAGroupFileEntry entry = getGroupFileEntry(aCycle);
	std::set<RNAGroupFileEntry>::const_iterator itEntry = mFiles.find(entry);
	assert(mFiles.end() != itEntry);

	CycleInfo cycle = aCycle;
	cycle.setChainAndOffset(' ', itEntry->offset());
	return cycle;
}

RNAGroupModel::cycle_pair RNAGroupModel::getGenericCyclePair(
	const CycleInfo& aCycle1,
	const CycleInfo& aCycle2) const
{
	CycleInfo cycle1 = getGenericCycle(aCycle1);
	CycleInfo cycle2 = getGenericCycle(aCycle2);

	cycle_pair cyclePair;
	if(cycle1 < cycle2)
	{
		cyclePair = cycle_pair(cycle1, cycle2);
	}else
	{
		cyclePair = cycle_pair(cycle2, cycle1);
	}
	return cyclePair;
}

InteractionInfo RNAGroupModel::getGenericInteraction(const InteractionInfo& aInteraction) const
{
	RNAGroupFileEntry entry = getGroupFileEntry(aInteraction);
	std::set<RNAGroupFileEntry>::const_iterator itEntry = mFiles.find(entry);
	assert(mFiles.end() != itEntry);

	InteractionInfo inter = aInteraction;
	inter.setChainAndOffset(' ', itEntry->offset());
	return inter;
}

// Takes a generic cycle
bool RNAGroupModel::internalIsNew(const CycleInfo& aCycle) const
{
	std::set<CycleInfo>::const_iterator it = mCycles.find(aCycle);
	return it == mCycles.end();
}

bool RNAGroupModel::isNew(
	const CycleInfo& aCycle1,
	const CycleInfo& aCycle2) const
{
	cycle_pair cyclePair = getGenericCyclePair(aCycle1, aCycle2);
	return internalIsNew(cyclePair);
}

bool RNAGroupModel::internalIsNew(const cycle_pair& aCyclePair) const
{
	std::set<cycle_pair>::const_iterator it = mCyclePairs.find(aCyclePair);
	return it == mCyclePairs.end();
}

void RNAGroupModel::addCyclePair(
	const CycleInfo& aCycle1,
	const CycleInfo& aCycle2)
{
	cycle_pair cyclePair = getGenericCyclePair(aCycle1, aCycle2);
	assert(internalIsNew(cyclePair));
	mCyclePairs.insert(cyclePair);
}

bool RNAGroupModel::isNew(const interaction_cycle_pair& aInteraction) const
{
	interaction_cycle_pair inter(
		getGenericInteraction(aInteraction.first),
		getGenericCycle(aInteraction.second),
		getGenericCycle(aInteraction.third));
	return internalIsNew(inter);
}

void RNAGroupModel::addInteraction(const interaction_cycle_pair& aInter)
{
	interaction_cycle_pair inter(
			getGenericInteraction(aInter.first),
			getGenericCycle(aInter.second),
			getGenericCycle(aInter.third));
	assert(internalIsNew(inter));
	mInteractionCyclePairs.insert(inter);
}

bool RNAGroupModel::internalIsNew(const interaction_cycle_pair& aInteraction) const
{
	std::set<interaction_cycle_pair>::const_iterator it;
	it = mInteractionCyclePairs.find(aInteraction);
	return it == mInteractionCyclePairs.end();
}

}; // namespace annotate

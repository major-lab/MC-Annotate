/*
 * RNAGroupModel.h
 *
 *  Created on: Apr 22, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_RNAGroupModel_H_
#define _annotate_RNAGroupModel_H_

#include "CycleInfo.h"
#include "InteractionInfo.h"
#include <list>
#include <set>

#include "RNAGroupFile.h"

namespace annotate {

class RNAGroupModel
{
public:
	typedef tuple3<InteractionInfo, CycleInfo, CycleInfo> interaction_cycle_pair;

	RNAGroupModel() {}
	RNAGroupModel(const std::list<RNAGroupFileEntry>& aFiles);

	bool isNew(const CycleInfo& aCycle) const;
	void addCycle(const CycleInfo& aCycle);

	bool isNew(const CycleInfo& aCycle1, const CycleInfo& aCycle2) const;
	void addCyclePair(const CycleInfo& aCycle1, const CycleInfo& aCycle2);

	// Interaction and involved pair of cycles
	bool isNew(const tuple3<InteractionInfo, CycleInfo, CycleInfo>& aInteraction) const;
	void addInteraction(const interaction_cycle_pair& aInteraction);
private:
	typedef std::pair<CycleInfo, CycleInfo> cycle_pair;
	std::set<CycleInfo> mCycles;
	std::set<RNAGroupFileEntry> mFiles;

	std::set<std::pair<CycleInfo, CycleInfo> > mCyclePairs;
	std::set<interaction_cycle_pair> mInteractionCyclePairs;

	RNAGroupFileEntry getGroupFileEntry(const CycleInfo& aCycle) const;
	RNAGroupFileEntry getGroupFileEntry(const InteractionInfo& aInter) const;
	CycleInfo getGenericCycle(const CycleInfo& aCycle) const;
	bool internalIsNew(const CycleInfo& aCycle) const;

	cycle_pair getGenericCyclePair(
		const CycleInfo& aCycle1,
		const CycleInfo& aCycle2) const;
	InteractionInfo getGenericInteraction(const InteractionInfo& aInteraction) const;

	bool internalIsNew(const cycle_pair& aCyclePair) const;
	bool internalIsNew(const interaction_cycle_pair& aInteraction) const;
};

}; // namespace annotate

#endif /* _annotate_RNAGroupModel_H_ */

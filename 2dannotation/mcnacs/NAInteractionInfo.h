/*
 * NAInteractionInfo.h
 *
 *  Created on: Nov 12, 2009
 *      Author: blanchmf
 */

#ifndef NAInteractionInfo_H_
#define NAInteractionInfo_H_

#include "InteractionInfo.h"
#include "CycleInfo.h"
#include <set>

class NAInteractionInfo : public annotate::InteractionInfo
{
public:
	NAInteractionInfo(const annotate::InteractionInfo& aInteraction);

	// ACCESS ------------------------------------------------------------------
	void addFivePrimeCycle(const annotate::CycleInfo& aCycle);
	void addThreePrimeCycle(const annotate::CycleInfo& aCycle);

	const std::set<annotate::CycleInfo>& fivePrimeCycles() const {return mFivePrimeCycles;}
	const std::set<annotate::CycleInfo>& threePrimeCycles() const {return mThreePrimeCycles;}

	std::string toString() const;

private:
	std::set<annotate::CycleInfo> mFivePrimeCycles;
	std::set<annotate::CycleInfo> mThreePrimeCycles;
};


#endif /* NAInteractionInfo_H_ */

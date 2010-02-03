#include "Cycle.h"
#include "AnnotateModel.h"

#include <cassert>
#include <sstream>
#include <memory>


namespace annotate
{
	Cycle::Cycle(const interactions_list& aInteractions)
	{
		assert(0 < aInteractions.size());

		mInteractions = aInteractions;

		mResidues = getOrderedResidues(mInteractions);

		updateProfile();

		assert(checkIntegrity());
	}

	Cycle::~Cycle()
	{
		mResidues.clear();
		mInteractions.clear();
	}

	std::list<mccore::ResId> Cycle::getOrderedResidues(
		const interactions_list& aInteractions) const
	{
		assert(0 < aInteractions.size());
		std::list<mccore::ResId> resids;
		mccore::ResId firstResId = aInteractions.front().fResId;
		resids.push_back(firstResId);

		interactions_list::const_iterator it;
		for(it = aInteractions.begin(); it != aInteractions.end(); ++ it)
		{
			assert(it->fResId == resids.back());
			if(it->rResId != firstResId)
			{
				resids.push_back(it->rResId);
			}
		}
		return resids;
	}

	void Cycle::updateProfile()
	{
		// Compute the profile
		int i = 0;
		// Adds the structure
		std::list<mccore::ResId>::const_iterator itPrev = mResidues.end();
		std::list<mccore::ResId>::const_iterator it;
		for(it = mResidues.begin(); it != mResidues.end(); ++ it)
		{
			if(it == mResidues.begin())
			{
				i = 1;
			}
			else if(!areResiduesLinked(*itPrev, *it))
			{
				mProfile.push_back(i);
				i = 1;
			}
			else
			{
				++ i;
			}
			itPrev = it;
		}
		mProfile.push_back(i);
	}

	bool Cycle::checkIntegrity() const
	{
		bool bIntegrity = true;
		interactions_list::const_iterator it;
		interactions_list::const_iterator itPrev;
		for(it = mInteractions.begin(); it != mInteractions.end() && bIntegrity; ++ it)
		{
			if(it != mInteractions.begin())
			{
				bIntegrity = (itPrev->rResId == it->fResId);
			}
			itPrev = it;
		}
		return bIntegrity;
	}

	bool Cycle::areResiduesLinked(
		const mccore::ResId& aRes1,
		const mccore::ResId& aRes2) const
	{
		bool bLinked = false;
		interactions_list::const_iterator it;
		for(it = mInteractions.begin(); it != mInteractions.end() && !bLinked; ++ it)
		{

			if(it->type() == BaseInteraction::eLINK
				&& ((aRes1 == it->fResId && aRes2 == it->rResId)
				|| (aRes1 == it->rResId && aRes2 == it->fResId)))
			{
				bLinked = true;
			}
		}
		return bLinked;
	}

	bool Cycle::areResiduesPaired(
		const mccore::ResId& aRes1,
		const mccore::ResId& aRes2) const
	{
		bool bPaired = false;
		interactions_list::const_iterator it;
		for(it = mInteractions.begin(); it != mInteractions.end() && !bPaired; ++ it)
		{
			if(it->type() == BaseInteraction::ePAIR
				&& ((aRes1 == it->fResId && aRes2 == it->rResId)
				|| (aRes1 == it->rResId && aRes2 == it->fResId)))
			{
				bPaired = true;
			}
		}
		return bPaired;
	}

	bool Cycle::areResiduesStacked(
		const mccore::ResId& aRes1,
		const mccore::ResId& aRes2) const
	{
		bool bStacked = false;
		interactions_list::const_iterator it;
		for(it = mInteractions.begin(); it != mInteractions.end() && !bStacked; ++ it)
		{
			if(it->type() == BaseInteraction::eSTACK
				&& ((aRes1 == it->fResId && aRes2 == it->rResId)
				|| (aRes1 == it->rResId && aRes2 == it->fResId)))
			{
				bStacked = true;
			}
		}
		return bStacked;
	}

	bool Cycle::operator <(const Cycle& aCycle) const
	{
		bool bSmaller = std::lexicographical_compare(
			mResidues.begin(),
			mResidues.end(),
			aCycle.mResidues.begin(),
			aCycle.mResidues.end());
		return bSmaller;
	}

	bool Cycle::shareInteractions(const Cycle& aCycle) const
	{
		bool bShare = false;
		std::set<BaseInteraction> thisInteractions;
		std::set<BaseInteraction> otherInteractions;

		thisInteractions.insert(mInteractions.begin(), mInteractions.end());
		otherInteractions.insert(aCycle.mInteractions.begin(), aCycle.mInteractions.end());

		std::set<BaseInteraction>::const_iterator first1 = thisInteractions.begin();
		std::set<BaseInteraction>::const_iterator last1 = thisInteractions.end();
		std::set<BaseInteraction>::const_iterator first2 = otherInteractions.begin();
		std::set<BaseInteraction>::const_iterator last2 = otherInteractions.end();
		while (first1!=last1 && first2!=last2 && !bShare)
  		{
			if (*first1<*first2)
			{
				 ++first1;
			}
			else if (*first2<*first1)
			{
				 ++first2;
			}
			else
			{
				bShare = true;
			}
		}
		return bShare;
	}

	bool Cycle::isSingleChain() const
	{
		bool bSingleChain = false;
		if(0 < mResidues.size())
		{
			bSingleChain = true;
			std::list<mccore::ResId>::const_iterator it = mResidues.begin();
			char chainId = it->getChainId();
			++ it;
			while(it != mResidues.end() && bSingleChain)
			{
				bSingleChain = (it->getChainId() == chainId);
				++ it;
			}
		}
		return bSingleChain;
	}

	bool Cycle::isClosed() const
	{
		bool bIsClosed = false;
		mccore::ResId frontId = mResidues.front();
		mccore::ResId backId = mResidues.back();

		if(	areResiduesPaired(frontId, backId)
			|| areResiduesStacked(frontId, backId)
			|| areResiduesLinked(frontId, backId))
		{
			bIsClosed = true;
		}
		return bIsClosed;
	}

	bool Cycle::isParallel() const
	{
		bool bIsParallel = false;

		if(2 == mProfile.size())
		{
			std::vector<std::vector<mccore::ResId> > strands = getStrands();
			assert(2 == strands.size());
			// TODO : Add support for 'triangle' connections
			assert(1 < strands[0].size());
			assert(1 < strands[1].size());

			if(	strands[0].front() < strands[0].back()
				&& strands[1].back() < strands[1].front())
			{
				bIsParallel = true;
			}else if( strands[0].back() < strands[0].front()
				&& strands[1].front() < strands[1].back())
			{
				bIsParallel = true;
			}
		}
		return bIsParallel;

	}

	Cycle::enType Cycle::getType() const
	{
		Cycle::enType eType;
		assert(0 < mProfile.size());
		if(mResidues.size() == 2 || !isClosed())
		{
			eType = eLOOSE;
		}else if(1 == mProfile.size())
		{
			eType = eLOOP;
		}else if(2 == mProfile.size())
		{
			std::vector<std::vector<mccore::ResId> > strands = getStrands();
			assert(2 == strands.size());
			if(1 == strands[0].size() || 1 == strands[1].size())
			{
				eType = e2STRANDS_TRIANGLE;
			}
			else if(isParallel())
			{
				eType = e2STRANDS_PARALLEL;
			}
			else
			{
				eType = e2STRANDS_ANTIPARALLEL;
			}
		}else
		{
			// Multibranch
			eType = eMULTIBRANCH;
		}
		return eType;
	}

	std::vector<std::vector<mccore::ResId> > Cycle::getStrands() const
	{
		std::vector<std::vector<mccore::ResId> > strands;
		strands.resize(mProfile.size());
		std::list<mccore::ResId>::const_iterator itRes = mResidues.begin();
		std::vector<unsigned int>::const_iterator itProf;
		for(unsigned int iStrand = 0; iStrand < mProfile.size(); ++ iStrand)
		{
			std::vector<mccore::ResId> strand;
			strand.resize(mProfile[iStrand]);
			for(unsigned int iRes = 0; iRes < mProfile[iStrand]; ++ iRes)
			{
				strand[iRes] = *itRes;
				++ itRes;
			}
			strands[iStrand] = strand;
		}
		return strands;
	}

	std::set<BaseInteraction> Cycle::getBaseInteractions() const
	{
		std::set<BaseInteraction> inters;
		std::list<BaseInteraction>::const_iterator it;
		for(it = mInteractions.begin();	it != mInteractions.end(); ++ it)
		{
			inters.insert(*it);
		}

		return inters;
	}
}

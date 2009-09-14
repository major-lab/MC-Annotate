#include "Cycle.h"
#include "AnnotateModel.h"
#include "mccore/AbstractModel.h"
#include "mccore/Pdbstream.h"

#include <cassert>
#include <sstream>
#include <memory>


namespace annotate
{
	Cycle::Cycle(const interactions_set& aInteractions)
	{
		assert(0 < aInteractions.size());

		setInteractions(aInteractions);

		mResidues = getOrderedResidues(mInteractions);

		updateProfile();
	}

	Cycle::~Cycle()
	{
		mResidues.clear();
		mInteractions.clear();
	}

	std::list<mccore::ResId> Cycle::getOrderedResidues(
		const interactions_set& aInteractions) const
	{
		std::list<mccore::ResId> resids;
		interactions_set interactions;
		interactions.insert(aInteractions.begin(), aInteractions.end());

		interactions_set::const_iterator itInter = interactions.begin();
		if(itInter != interactions.end())
		{
			resids.push_back(itInter->fResId);
			resids.push_back(itInter->rResId);
			interactions.erase(itInter);
		}

		unsigned int iSize = interactions.size();
		while(0 < iSize)
		{
			interactions_set::const_iterator it;
			for(it = interactions.begin(); it != interactions.end(); ++ it)
			{
				if(it->fResId == resids.front())
				{
					resids.push_front(it->rResId);
					interactions.erase(it);
					break;
				}
				else if(it->rResId == resids.front())
				{
					resids.push_front(it->fResId);
					interactions.erase(it);
					break;
				}
				else if(it->fResId == resids.back())
				{
					resids.push_back(it->rResId);
					interactions.erase(it);
					break;
				}
				else if(it->rResId == resids.back())
				{
					resids.push_back(it->fResId);
					interactions.erase(it);
					break;
				}
			}
			assert(interactions.size() < iSize);
			iSize = interactions.size();
		}

		// Complete cycle, remove duplicates
		if(resids.front() == resids.back())
		{
			resids.pop_back();
		}

		// Reorder as necessary
		orderResidues(aInteractions, resids);

		return resids;
	}

	bool Cycle::isClosed(
		const interactions_set& aInteractions,
		std::list<mccore::ResId>& aResidues) const
	{
		mccore::ResId minResId = std::min(aResidues.front(), aResidues.back());
		mccore::ResId maxResId = std::max(aResidues.front(), aResidues.back());

		interactions_set::const_iterator it;
		for(it = aInteractions.begin(); it != aInteractions.end(); ++ it)
		{
			if((std::min(it->fResId, it->rResId) == minResId)
			&& (std::max(it->fResId, it->rResId) == maxResId))
			{
				break;
			}
		}

		return (it != aInteractions.end());
	}

	const Cycle::interactions_set Cycle::getInteractions() const
	{
		interactions_set interactions;
		interactions.insert(mInteractions.begin(), mInteractions.end());
		return interactions;
	}

	void Cycle::orderResidues(
		const interactions_set& aInteractions,
		std::list<mccore::ResId>& aResidues) const
	{
		BaseInteraction inter(0, aResidues.front(), 0, aResidues.back());

		if(isClosed(aInteractions, aResidues))
		{
			// Complete cycle, rotate until the first residue is the lowest
			std::list<mccore::ResId>::iterator itMin;
			itMin = std::min_element(aResidues.begin(), aResidues.end());
			std::rotate(aResidues.begin(), itMin, aResidues.end());

			// Order the residue so the cycle is "slowly growing"
			std::list<mccore::ResId>::const_iterator it = aResidues.begin();
			it ++;
			if(aResidues.back() < *it)
			{
				aResidues.reverse();
				aResidues.push_front(aResidues.back());
				aResidues.pop_back();
			}
		}
		else if(aResidues.back() < aResidues.front())
		{
			// Insure the first residue is the minimum
			aResidues.reverse();
		}
	}

	void Cycle::setInteractions(const interactions_set& aInteractions)
	{
		interactions_set::const_iterator it;
		for(it = aInteractions.begin(); it != aInteractions.end(); ++ it)
		{
			assert(it->type() != BaseInteraction::eUNKNOWN);
			mInteractions.insert(*it);
		}
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

	bool Cycle::areResiduesLinked(
		const mccore::ResId& aRes1,
		const mccore::ResId& aRes2) const
	{
		bool bLinked = false;
		interactions_set::const_iterator it;
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
		interactions_set::const_iterator it;
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
		interactions_set::const_iterator it;
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
		interactions_set interactions;
		interactions_set::const_iterator first1 = mInteractions.begin();
		interactions_set::const_iterator last1 = mInteractions.end();
		interactions_set::const_iterator first2 = aCycle.mInteractions.begin();
		interactions_set::const_iterator last2 = aCycle.mInteractions.end();
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

	Cycle::interactions_set Cycle::getBaseInteractions() const
	{
		interactions_set inters;
		interactions_set::const_iterator it;
		for(it = mInteractions.begin();	it != mInteractions.end(); ++ it)
		{
			BaseInteraction inter = *it;
			inter.type() = BaseInteraction::eUNKNOWN;
			inters.insert(inter);
		}

		return inters;
	}
}

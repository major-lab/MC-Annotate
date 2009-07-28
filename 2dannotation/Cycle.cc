#include "Cycle.h"
#include "mccore/AbstractModel.h"
#include "mccore/Pdbstream.h"

// TODO : Remove this ( DEBUGGING )
#include "mccore/Messagestream.h"

#include <sstream>

namespace annotate
{	
	Cycle::Cycle(const mccore::GraphModel& aModel, unsigned char aucRelationMask)
	: mName("")
	{
		mucRelationMask = aucRelationMask;
		mModel = aModel;
		
		update();
	}
	
	Cycle::Cycle(
		const mccore::GraphModel& aModel, 
		const std::set<mccore::ResId>& aResIds,
		unsigned char aucRelationMask)
	{
		// Create a cycle from a set of base interactions
		mucRelationMask = aucRelationMask;
		mModel.clear();
		
		std::set<mccore::ResId>::const_iterator itRes;
		for(itRes = aResIds.begin(); itRes != aResIds.end(); ++ itRes)
		{
			GraphModel::const_iterator itResidue = aModel.find(*itRes);
			mModel.insert(*itResidue);
		}		
		update();
	} 
	
	Cycle::~Cycle()
	{
		mInteractions.clear();
		mResidues.clear();
	}
		
	const mccore::GraphModel& Cycle::getModel() const
	{
		return mModel;
	}
	
	void Cycle::update()
	{
		mModel.annotate(mucRelationMask);
			
		// Get the interactions and residues
		mInteractionsAnnotation.update(mModel);
		AbstractModel::const_iterator it;
		std::list<const BaseInteraction*> interactions;
		for(it = mModel.begin(); it != mModel.end(); ++it)
		{
			AbstractModel::const_iterator itNext = it;
			itNext ++;
			if(itNext == mModel.end())
			{
				itNext = mModel.begin();
			}
			mccore::ResId refId = it->getResId();
			mccore::ResId resId = itNext->getResId();			
			
			// TODO: Assert that this is not NULL
			interactions = mInteractionsAnnotation.getInteractions(refId, resId);
			mInteractions.insert(interactions.begin(), interactions.end());
			mResidues.push_back(refId);
		}
		
		updateProfile();
	}
	
	void Cycle::updateProfile()
	{
		std::list<mccore::ResId>::iterator itMin;
		std::list<mccore::ResId>::iterator itNext;
		itMin = std::min_element(mResidues.begin(), mResidues.end());
		std::rotate(mResidues.begin(), itMin, mResidues.end());
		itMin = mResidues.begin();
		itNext = itMin;
		itNext ++;
		
		// Insure first relation is increasing
		if(mResidues.back() < (*itNext))
		{
			std::reverse(mResidues.begin(), mResidues.end());
			mResidues.push_front(mResidues.back());
			mResidues.pop_back();
		}
		
		// Now the cycle is in order, starting with the 5' base
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
			else if(!itPrev->areContiguous(*it))
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
		interactions_set_iterator first1 = mInteractions.begin();
		interactions_set_iterator last1 = mInteractions.end();
		interactions_set_iterator first2 = aCycle.mInteractions.begin();
		interactions_set_iterator last2 = aCycle.mInteractions.end();	
		while (first1!=last1 && first2!=last2 && !bShare)
  		{
  			const BaseInteraction* pLeft = *first1;
  			const BaseInteraction* pRight = *first2;
  			
			if (*pLeft<*pRight)
			{
				 ++first1;
			}
			else if (*pRight<*pLeft)
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
	
	std::string Cycle::getSequence() const
	{
		std::ostringstream oss;
		std::list<mccore::ResId>::const_iterator it;
		for(it = mResidues.begin(); it != mResidues.end(); ++ it)
		{
			mccore::GraphModel::const_iterator itRes = mModel.find(*it);
			oss << mccore::Pdbstream::stringifyResidueType (itRes->getType());
		}
		return oss.str();
	}
	
	std::set<BaseInteraction> Cycle::getBaseInteractions() const
	{
		std::set<BaseInteraction> interactions;
		interactions_set_iterator it;
		// Get the set of resId to resId interactions, without qualification
		for(it = mInteractions.begin(); it != mInteractions.end(); ++ it)
		{
			const BaseInteraction* pBaseInter = *it;
			BaseInteraction inter(
				pBaseInter->first, pBaseInter->fResId, 
				pBaseInter->second, pBaseInter->rResId);
			interactions.insert(inter);	
		}		
		return interactions;
	}
}
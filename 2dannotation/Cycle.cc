#include "Cycle.h"
#include "AnnotateModel.h"
#include "mccore/AbstractModel.h"
#include "mccore/Pdbstream.h"

#include <sstream>
#include <cassert>

namespace annotate
{	
	Cycle::Cycle(const mccore::GraphModel& aModel, unsigned char aucRelationMask)
	: mName("")
	{
		// Create from a cycle model
		mucRelationMask = aucRelationMask;
		mModel = aModel;
		
		mModel.annotate(mucRelationMask);

		// Verify what we have for residues
		std::set<mccore::ResId> residues;
		GraphModel::const_iterator it;
		for(it = mModel.begin(); it != mModel.end(); ++it)
		{
			residues.insert(it->getResId());
		}
		
		// Get the interactions forming the cycle
		AnnotationInteractions annInteractions;
		annInteractions.update(mModel);
		std::list<const BaseInteraction*> interactions;
		interactions = annInteractions.getInteractions(residues);		
		setInteractions(interactions);
		
		mResidues = getOrderedResidues(mInteractions);
		
		updateProfile();
	}
	
	Cycle::Cycle(
		const AnnotateModel& aModel, 
		const std::set<BaseInteraction>& aInteractions,
		unsigned char aucRelationMask)
	{
		assert(0 < aInteractions.size());
		
		// Create a cycle from a set of interactions from another model
		mucRelationMask = aucRelationMask;
		
		mInteractions = aInteractions;
		mResidues = getOrderedResidues(mInteractions);
		
		// Insert the residues in the model
		std::list<mccore::ResId>::const_iterator it;
		for(it = mResidues.begin(); it != mResidues.end(); ++it)
		{
			mccore::GraphModel::const_iterator itRes = aModel.find(*it);
			mModel.insert(*itRes);
		}
		
		updateProfile();
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
	
	mccore::ResId Cycle::getMinResId(const interactions_set& aInteractions) const
		throw(mccore::NoSuchElementException)
	{
		mccore::ResId minResId;
		if(0 < aInteractions.size())
		{
			interactions_set::const_iterator it = aInteractions.begin();
			minResId = std::min(it->fResId, it->rResId);
			for(; it != aInteractions.end(); ++ it)
			{
				if(it->fResId < minResId)
				{
					minResId = it->fResId;
				}
				if(it->rResId < minResId)
				{
					minResId = it->rResId;
				}
			}
		}
		else
		{
			std::string strMsg("Cycle::getMinResId - Interaction set empty");
			throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
		}
		return minResId;
	}
	
	std::list<mccore::ResId> Cycle::getOrderedResidues(
		const interactions_set& aInteractions) const
	{
		std::list<mccore::ResId> resids;
		interactions_set interactions = aInteractions;
		mccore::ResId resId = getMinResId(interactions);
		
		// Get all the residues ids from the interactions
		unsigned int iSize = interactions.size();
		while(0 < iSize)
		{
			interactions_set::const_iterator it;
			for(it = interactions.begin(); it != interactions.end(); ++ it)
			{
				if(resId == it->fResId)
				{
					resids.push_back(resId);
					resId = it->rResId;
					interactions.erase(it);
					break;					
				}
				else if(resId == it->rResId)
				{
					resids.push_back(resId);
					resId = it->fResId;
					interactions.erase(it);
					break;
				}
			}
			assert(interactions.size() < iSize);
			iSize = interactions.size();
		}
		
		// Reorder as necessary
		std::list<mccore::ResId>::const_iterator it = resids.begin();
		it ++;
		if(resids.back() < *it)
		{
			resids.reverse();
			resids.push_front(resids.back());
			resids.pop_back();
		}
		return resids;
	}
	
	void Cycle::setInteractions(const std::list<const BaseInteraction*>& aInteractions)
	{
		mInteractions.clear();
		std::list<const BaseInteraction*>::const_iterator it;
		for(it = aInteractions.begin(); it != aInteractions.end(); ++ it)
		{
			const BaseInteraction* pInter = *it;			
			BaseInteraction interact(
				pInter->first, 
				pInter->fResId, 
				pInter->second, 
				pInter->rResId);
			mInteractions.insert(interact);
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
	
	std::string Cycle::getSequence() const
		throw(mccore::NoSuchElementException)
	{
		std::ostringstream oss;
		std::list<mccore::ResId>::const_iterator it;
		for(it = mResidues.begin(); it != mResidues.end(); ++ it)
		{
			mccore::GraphModel::const_iterator itRes = mModel.find(*it);
			assert(itRes != mModel.end()); // Cycle is invalid
			oss << mccore::Pdbstream::stringifyResidueType (itRes->getType());
		}
		return oss.str();
	}
}
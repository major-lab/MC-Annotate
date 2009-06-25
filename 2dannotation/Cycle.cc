#include "Cycle.h"
#include "mccore/AbstractModel.h"

// TODO : Remove this ( DEBUGGING )
#include "mccore/Messagestream.h"

namespace annotate
{	
	Cycle::Cycle(const mccore::GraphModel& aModel)
	{
		mModel = aModel;
	}
	
	Cycle::~Cycle()
	{
	}
		
	const mccore::GraphModel& Cycle::getModel() const
	{
		return mModel;
	}
	
	void Cycle::order()
	{
		mModel.annotate();
		
		// Get the interactions and residues
		mInteractionsAnnotation.update(mModel);
		AbstractModel::const_iterator it;
		for(it = mModel.begin(); it != mModel.end(); ++it)
		{
			const BaseInteraction* pInteraction;
			AbstractModel::const_iterator itNext = it;
			itNext ++;
			if(itNext == mModel.end())
			{
				itNext = mModel.begin();
			}
			mccore::ResId resId = it->getResId();
			mccore::ResId refId = itNext->getResId();			
			
			// TODO: Assert that this is not NULL
			pInteraction = mInteractionsAnnotation.getInteraction(std::min(resId, refId), std::max(resId, refId));
			if(NULL == pInteraction)
			{
				gOut (0) << "No interaction found between " << resId << " and " << refId << std::endl;
			}
			mInteractions.push_back(pInteraction);
		}
		
		mTopology = computeTopology();
	}
	
	const std::list<unsigned int>& Cycle::getTopology() const
	{
		return mTopology.getStrands();
	}
	
	Cycle::Topology Cycle::computeTopology() const
	{
		Topology topo;
		interactions_list_iterator it;
		for(it = mInteractions.begin(); it != mInteractions.end(); ++it)
		{
			const BaseInteraction* pInteraction = *it;
			const BasePair* pPair = dynamic_cast<const BasePair*>(pInteraction);
			if(NULL != pPair)
			{
				Topology candidate(it, false);
				completeTopology(candidate);
				if(topo.empty() || topo < candidate)
				{
					topo = candidate;					
				}
				
				Topology reverseCandidate(it, true);
				completeReverseTopology(reverseCandidate);
				if(topo < reverseCandidate)
				{
					topo = reverseCandidate;					
				}
			}
		}
		return topo;
	}
	
	void Cycle::completeTopology(Topology& aCandidate) const
	{
		Topology::const_interaction_iterator it = aCandidate.getFirstIterator();
		Topology::const_interaction_iterator itFirst = it;
		do
		{
			unsigned int iCount = advanceTopoIterator(it);
			aCandidate.addStrandLength(iCount);
		} while(it != itFirst);
	}
	
	void Cycle::completeReverseTopology(Topology& aCandidate) const
	{
		Topology::const_interaction_iterator it = aCandidate.getFirstIterator();
		Topology::const_interaction_iterator itFirst = it;
		do
		{
			unsigned int iCount = rewindTopoIterator(it);
			aCandidate.addStrandLength(iCount);
		} while(it != itFirst);
	}
	
	unsigned int Cycle::advanceTopoIterator(interactions_list_iterator& it) const
	{
		unsigned int iCount = 0;
		const BaseInteraction* pPair = NULL;
		while(NULL == pPair)
		{
			it ++;
			iCount ++;
			if(it == mInteractions.end())
			{
				it = mInteractions.begin();
			}
			pPair = dynamic_cast<const BasePair*>(*it);
		} 
		return iCount;
	}
	
	unsigned int Cycle::rewindTopoIterator(interactions_list_iterator& it) const
	{
		unsigned int iCount = 0;
		const BaseInteraction* pPair = NULL;
		while(NULL == pPair)
		{
			if(it == mInteractions.begin())
			{
				it = mInteractions.end();
			}
			it --;
			
			pPair = dynamic_cast<const BasePair*>(*it);
			if(NULL == pPair)
			{
				iCount ++;
			}
		} 
		return iCount;
	}
	
	bool Cycle::Topology::operator <(const Topology& aRight) const
	{
		bool bSmaller = false;
		std::list<unsigned int>::const_iterator it1 = mStrands.begin();
		std::list<unsigned int>::const_iterator it2 = aRight.mStrands.begin();
		if(mStrands.size() == aRight.mStrands.size())
		{
			for(; it1 != mStrands.end() && it2 != aRight.mStrands.end() && !bSmaller; ++it1, ++it2)
			{
				bSmaller = (*it1 < *it2);
			}
		}
		return bSmaller;
	}
}
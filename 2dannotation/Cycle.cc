#include "Cycle.h"
#include "mccore/AbstractModel.h"
#include "mccore/Pdbstream.h"

// TODO : Remove this ( DEBUGGING )
#include "mccore/Messagestream.h"

namespace annotate
{	
	Cycle::Cycle(const mccore::GraphModel& aModel, unsigned char aucRelationMask)
	: mName("")
	{
		mucRelationMask = aucRelationMask;
		mModel = aModel;
		
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
			
			if(0 == interactions.size())
			{
				gOut (0) << "No interaction found between " << resId << " and " << refId << std::endl;
			}
			mInteractions.insert(interactions.begin(), interactions.end());
			mResidues.push_back(refId);
		}
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
}
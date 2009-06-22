#include "Cycle.h"
#include "AnnotationInteractions.h"
#include "mccore/AbstractModel.h"

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
		// Get the interactions and residues
		AnnotationInteractions annot;
		annot.update(mModel);
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
			// TODO: Assert that this is not NULL
			mccore::ResId resId = it->getResId();
			mccore::ResId refId = itNext->getResId();
			pInteraction = annot.getInteraction(resId, refId);
			mInteractions.push_back(pInteraction);
			mResidues.push_back(&(*it));
		}
	}
}
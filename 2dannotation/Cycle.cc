#include "Cycle.h"
#include "mccore/AbstractModel.h"
#include "mccore/Pdbstream.h"

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
		mInteractions.clear();
		mResidues.clear();
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
		interactions_list interactions;
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
			mInteractions.insert(mInteractions.end(), interactions.begin(), interactions.end());
			mResidues.push_back(&(*it));
		}
	}
	
	std::list<std::string> Cycle::getInteractionLabels() const
	{
		std::list<std::string> labels;
		std::string currentLabel;
		std::list<const BaseInteraction*>::const_iterator it = mInteractions.begin();
		std::list<const mccore::Residue*>::const_iterator itRes = mResidues.begin();
		const BaseInteraction* prev = *it;
		for(; it != mInteractions.end() && itRes != mResidues.end(); ++ it, ++itRes)
		{
			if(!prev->sameResidues(*(*it)))
			{
				// const mccore::Residue* prevResidue = &(*mModel.find((*prev).fResId));
				//std::string strResidue = mccore::Pdbstream::stringifyResidueType (prevResidue->getType());
				//currentLabel = strResidue + currentLabel;
				labels.push_back(currentLabel);
				currentLabel = "";
			}
			prev = *it;
			
			if(NULL != dynamic_cast<const BasePair*>(*it))
			{
				currentLabel += "P";
			}else if(NULL != dynamic_cast<const BaseLink*>(*it))
			{
				currentLabel = "L" + currentLabel;
			}else if(NULL != dynamic_cast<const BaseStack*>(*it))
			{
				currentLabel = "S" + currentLabel;
			}			
		}
		if(!currentLabel.empty())
		{
			// const mccore::Residue* prevResidue = &(*mModel.find((*prev).fResId));
			// std::string strResidue = mccore::Pdbstream::stringifyResidueType (prevResidue->getType());
			// currentLabel = strResidue + currentLabel;
			labels.push_back(currentLabel);
		}
		return labels;
	}
/*	
	std::pair<bool, unsigned int> Cycle::epitomeTranformation() const
	{
		bool bReverse;
		unsigned int uiRotation;
	}
	
	unsigned int Cycle::epitomeTranformation(std::list<std::string> interactions) const
	{
		std::list<unsigned int> interactScores;
		std::list<std::string>::const_iterator it;
		for(it = interactions.begin(); it != interactions.end(); ++it)
		{
			int score = 0;
			if(string::npos != it->find('L'))
			{
				iScore += 4;
			}
			if(string::npos != it->find('P'))
			{
				iScore += 2;
			}
			if(string::npos != it->find('S'))
			{
				iScore += 1;
			}
			interactScores.push_back(iScore);
		}
	}
	*/
}
#include "AnnotationTertiaryCycles.h"
#include "AnnotationCycles.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotationTertiaryStacks.h"
#include "AnnotationInteractions.h"
#include "AnnotateModel.h"
#include "AlgorithmExtra.h"
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationTertiaryCycles::mstrAnnotationName = "Tertiary Cycles";
	
	// Methods	
	AnnotationTertiaryCycles::AnnotationTertiaryCycles()
	{
		addRequirement<AnnotationCycles>();
		addRequirement<AnnotationTertiaryPairs>();
		addRequirement<AnnotationTertiaryStacks>();
	}
	
	AnnotationTertiaryCycles::~AnnotationTertiaryCycles()
	{
		clear();
	}
	
	void AnnotationTertiaryCycles::clear()
	{
		mCycles.clear();
	}
		
	void AnnotationTertiaryCycles::update(const AnnotateModel& aModel)
	{
		clear();
		
		const AnnotationCycles* pAnnCycles = 
			aModel.getAnnotation<AnnotationCycles>();
		const AnnotationTertiaryPairs* pAnnTertiaryPairs = 
			aModel.getAnnotation<AnnotationTertiaryPairs>();
		const AnnotationTertiaryStacks* pAnnTertiaryStacks = 
			aModel.getAnnotation<AnnotationTertiaryStacks>();
		
		if(NULL != pAnnCycles && NULL != pAnnTertiaryPairs && NULL != pAnnTertiaryStacks)
		{
			std::set<BaseInteraction> cyclePairs;
			mCycles = pAnnCycles->getCycles();
			bool bTertiary = false;
			std::list<Cycle>::iterator it = mCycles.begin();
			int i = 0;
			while(it != mCycles.end())
			{
				bTertiary = false;
				if(it->isSingleChain())
				{
					// Check if the cycle is tertiary
					getPairs(*it, cyclePairs);
				
					bTertiary = isTertiary(
						cyclePairs, 
						pAnnTertiaryPairs->getPairs());
					if(!bTertiary)
					{
						bTertiary = isTertiary(
							cyclePairs, 
							pAnnTertiaryStacks->getStacks());
					}
					cyclePairs.clear();
				}
				
				if(bTertiary)
				{
					++it;
				}
				else
				{
					it = mCycles.erase(it);
				}
				++ i;
				
			}
		}
	}
	
	std::string AnnotationTertiaryCycles::output() const
	{
		std::ostringstream oss;
		std::list<Cycle>::const_iterator itCycle;
		for (itCycle = mCycles.begin(); mCycles.end() != itCycle; ++itCycle)
	    {
	    	oss << itCycle->name() << std::endl;
	    }	  	
		return oss.str();
	}
	
	bool AnnotationTertiaryCycles::isTertiary(
		const std::set<BaseInteraction>& aCyclePairs, 
		const std::set<BasePair>& a3DPairs) const
	{
		bool bTertiary = false;
		std::set<BaseInteraction>::const_iterator first1 = aCyclePairs.begin();
		std::set<BaseInteraction>::const_iterator last1 = aCyclePairs.end();
		std::set<BasePair>::const_iterator first2 = a3DPairs.begin();
		std::set<BasePair>::const_iterator last2 = a3DPairs.end();
		while (first1!=last1 && first2!=last2)
  		{
  			const BaseInteraction* pLeft = &(*first1);
  			const BaseInteraction* pRight = &(*first2);
  			
			if (*pLeft<*pRight) ++first1;
			else if (*pRight<*pLeft) ++first2;
			else 
			{
				bTertiary = true;
				break; 
			}
		}				
		return bTertiary;		
	}
	
	bool AnnotationTertiaryCycles::isTertiary(
		const std::set<BaseInteraction>& aCyclePairs, 
		const std::set<BaseStack>& a3DStacks) const
	{
		bool bTertiary = false;
		std::set<BaseInteraction>::const_iterator first1 = aCyclePairs.begin();
		std::set<BaseInteraction>::const_iterator last1 = aCyclePairs.end();
		std::set<BaseStack>::const_iterator first2 = a3DStacks.begin();
		std::set<BaseStack>::const_iterator last2 = a3DStacks.end();
		while (first1!=last1 && first2!=last2)
  		{
  			const BaseInteraction* pLeft = &(*first1);
  			const BaseInteraction* pRight = &(*first2);
  			
			if (*pLeft<*pRight) ++first1;
			else if (*pRight<*pLeft) ++first2;
			else 
			{
				bTertiary = true;
				break; 
			}
		}				
		return bTertiary;		
	}
	
	void AnnotationTertiaryCycles::getPairs(
		const Cycle& aCycle, 
		std::set<BaseInteraction>& aPairs) const
	{		
    	AbstractModel::const_iterator it;
		for (it = aCycle.getModel().begin(); aCycle.getModel().end() != it; ++it)
		{
			mccore::ResId resId1 = it->getResId();
			mccore::ResId resId2;
			AbstractModel::const_iterator itNext = it;
			itNext ++;
			if(itNext == aCycle.getModel().end())
			{
				resId2 = aCycle.getModel().front().getResId();
			}
			else
			{
				resId2 = itNext->getResId();
			}
			if(resId2 < resId1)
			{
				std::swap(resId1, resId2);
			}
			BaseInteraction interaction(0, resId1, 0, resId2);
			aPairs.insert(interaction);
		}
	}
}
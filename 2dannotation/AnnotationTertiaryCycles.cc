#include "AnnotationTertiaryCycles.h"
#include "AnnotationCycles.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotateModel.h"
#include "AlgorithmExtra.h"
#include <sstream>

namespace annotate
{	
	AnnotationTertiaryCycles::AnnotationTertiaryCycles()
	{
		addRequirement(AnnotationCycles().provides());
		addRequirement(AnnotationTertiaryPairs().provides());
	}
	
	AnnotationTertiaryCycles::~AnnotationTertiaryCycles()
	{
		clear();
	}
	
	void AnnotationTertiaryCycles::clear()
	{
		mCycles.clear();
	}
	
	const std::string AnnotationTertiaryCycles::provides() const
	{
		std::string strAnnotationName = "Tertiary Cycles";
		return strAnnotationName;
	}
		
	void AnnotationTertiaryCycles::update(const AnnotateModel& aModel)
	{
		clear();
		
		const AnnotationCycles* pAnnCycles = 
			aModel.getAnnotation<AnnotationCycles>(AnnotationCycles().provides());
		const AnnotationTertiaryPairs* pAnnTertiaryPairs = 
			aModel.getAnnotation<AnnotationTertiaryPairs>(AnnotationTertiaryPairs().provides());
		
		if(NULL != pAnnCycles && NULL != pAnnTertiaryPairs)
		{
			std::set<BaseInteraction> cyclePairs;
			mCycles = pAnnCycles->getCycles();
			bool bTertiary = false;
			std::vector<Cycle>::iterator it = mCycles.begin();
			int i = 0;
			while(it != mCycles.end())
			{
				// Check if the cycle is tertiary
				getPairs(*it, cyclePairs);
				bTertiary = isTertiary(
					cyclePairs, 
					pAnnTertiaryPairs->getPairs());
				cyclePairs.clear();
				
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
	
	void AnnotationTertiaryCycles::outputCycle(
		std::ostringstream& oss, 
		const Cycle& aCycle) const
	{
		AbstractModel::const_iterator it;
    	oss << "[";
		for (it = aCycle.getModel().begin(); aCycle.getModel().end() != it; ++it)
		{
			oss << " " << it->getResId () << *it->getType ();
		}
		oss << " ] " << aCycle.getModel().size() << endl;
	}
	
	std::string AnnotationTertiaryCycles::output() const
	{
		int i = 0;
		std::ostringstream oss;
		std::vector<Cycle>::const_iterator itCycle;
		for (itCycle = mCycles.begin(); mCycles.end() != itCycle; ++itCycle)
	    {
	    	oss << "Cycle " << i << " : ";
	    	outputCycle(oss, *itCycle);
	    	++ i;
	    }	  	
		return oss.str();		
	}
	
	const std::vector< Cycle >& AnnotationTertiaryCycles::getCycles() const
	{
		return mCycles;
	}
	
	std::vector< Cycle >& AnnotationTertiaryCycles::getCycles()
	{
		return mCycles;
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
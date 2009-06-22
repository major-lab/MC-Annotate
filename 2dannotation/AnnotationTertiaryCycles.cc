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
			resid_pair_set cyclePairs;
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
	
	bool AnnotationTertiaryCycles::compare_less(
		const BasePair& aLeft, 
		const resid_pair& aRight) const
	{
		mccore::ResId leftMin = std::min(aLeft.fResId, aLeft.rResId);
		mccore::ResId rightMin = std::min(aRight.first, aRight.second);
		return ( leftMin < rightMin 
			|| (leftMin == rightMin 
				&& (std::max(aLeft.fResId, aLeft.rResId) < std::max(aRight.first, aRight.second))));
	}
	
	bool AnnotationTertiaryCycles::compare_less(
		const resid_pair& aLeft, 
		const BasePair& aRight) const
	{
		mccore::ResId rightMin = std::min(aRight.fResId, aRight.rResId);
		mccore::ResId leftMin = std::min(aLeft.first, aLeft.second);
		return ( leftMin < rightMin 
			|| (leftMin == rightMin 
				&& (std::max(aLeft.first, aLeft.second) < std::max(aRight.fResId, aRight.rResId))));
	}
	
	bool AnnotationTertiaryCycles::isTertiary(
		const resid_pair_set& aCyclePairs, 
		const std::vector<BasePair>& a3DPairs) const
	{		
		bool bTertiary = false;
		
		resid_pair_set::const_iterator first1 = aCyclePairs.begin();
		resid_pair_set::const_iterator last1 = aCyclePairs.end();
		std::vector<BasePair>::const_iterator first2 = a3DPairs.begin();
		std::vector<BasePair>::const_iterator last2 = a3DPairs.end();
			
		while (first1!=last1 && first2!=last2)
		{
		    if (compare_less(*first1, *first2))
		    {
		    	 ++first1;
		    }
		    else if (compare_less(*first2, *first1))
		    {
		    	 ++first2;
		    }
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
		resid_pair_set& aPairs) const
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
			if(resId1 < resId2)
			{
				std::swap(resId1, resId2);
			}
			std::pair<mccore::ResId, mccore::ResId> pair(resId1, resId2);
			aPairs.insert(pair);
		}
	}
}
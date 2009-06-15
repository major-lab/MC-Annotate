#include "AnnotationTertiaryCycles.h"
#include "AnnotationCycles.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotateModel.h"
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
			mCycles = pAnnCycles->getCycles();
			bool bTertiary = false;
			std::vector<Cycle>::iterator it = mCycles.begin();
			int i = 0;
			while(it != mCycles.end())
			{
				// Check if the cycle is tertiary
				bTertiary = isTertiary(
					getPairs(*it), 
					getPairs(*pAnnTertiaryPairs));
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
		for (it = aCycle.begin(); aCycle.end() != it; ++it)
		{
			oss << " " << it->getResId () << *it->getType ();
		}
		oss << " ] " << aCycle.size() << endl;
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
	
	bool AnnotationTertiaryCycles::isTertiary(
		const std::pair<mccore::ResId, mccore::ResId>& aCyclePair,
		const resid_pair_vector& a3DPairs) const
	{
		bool bTertiary = false;
		resid_pair_vector::const_iterator it = a3DPairs.begin();
		it = std::find(a3DPairs.begin(), a3DPairs.end(), aCyclePair);
		bTertiary = (it != a3DPairs.end());
		return bTertiary;
	}
	
	bool AnnotationTertiaryCycles::isTertiary(
		const resid_pair_vector& aCyclePairs, 
		const resid_pair_vector& a3DPairs) const
	{		
		bool bTertiary = false;
		resid_pair_vector::const_iterator itCP;
		for(itCP = aCyclePairs.begin(); 
			itCP != aCyclePairs.end() && !bTertiary; 
			++itCP)
		{
			bTertiary = isTertiary(*itCP, a3DPairs);
		}
		return bTertiary;		
	}
	
	AnnotationTertiaryCycles::resid_pair_vector
	AnnotationTertiaryCycles::getPairs(const Cycle& aCycle) const
	{
		resid_pair_vector pairs;    		
    	AbstractModel::const_iterator it;
		for (it = aCycle.begin(); aCycle.end() != it; ++it)
		{
			mccore::ResId resId1 = it->getResId();
			mccore::ResId resId2;
			AbstractModel::const_iterator itNext = it;
			itNext ++;
			if(itNext == aCycle.end())
			{
				resId2 = aCycle.front().getResId();
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
			
		
			pairs.push_back(pair);
		}
		return pairs;
	}
	
	AnnotationTertiaryCycles::resid_pair_vector 
	AnnotationTertiaryCycles::getPairs(
		const AnnotationTertiaryPairs& aTPAnnotation) const
	{
		resid_pair_vector pairs;
		std::vector<BasePair>::const_iterator it;
		for(
			it = aTPAnnotation.getPairs().begin(); 
			it != aTPAnnotation.getPairs().end(); 
			++ it)
		{
			ResId resId1 = it->fResId;
			ResId resId2 = it->rResId;
			
			if(resId1 < resId2)
			{
				std::swap(resId1, resId2);
			}
			std::pair<mccore::ResId, mccore::ResId> pair(resId1, resId2);
			pairs.push_back(pair);			
		}
		return pairs;
	}
}
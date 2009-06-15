#include "AnnotationTertiaryPairs.h"
#include "AnnotateModel.h"
#include "AnnotationStems.h"
#include "AnnotationLoops.h"
#include <sstream>

namespace annotate
{
	AnnotationTertiaryPairs::AnnotationTertiaryPairs()
	{
		// Requires 
		addRequirement(AnnotationStems().provides());
		addRequirement(AnnotationLoops().provides());
	}
	
	AnnotationTertiaryPairs::~AnnotationTertiaryPairs()
	{
		clear();
	}
	
	void AnnotationTertiaryPairs::clear()
	{
		mPairs.clear();
	}
	
	const std::string AnnotationTertiaryPairs::provides() const
	{
		std::string strAnnotationName = "Tertiary Pairs";
		return strAnnotationName;
	}
		
	const std::vector< BasePair >& AnnotationTertiaryPairs::getPairs() const
	{
		return mPairs;
	}
	
	void AnnotationTertiaryPairs::update(const AnnotateModel& aModel)
	{
		struct_association_map associationMap;
		
		// Add the stem associations
		addStemsAssociations(aModel, associationMap);
		
		// Add the loop associations
		addLoopsAssociations(aModel, associationMap);
		
		std::vector<BasePair>::const_iterator it = aModel.getBasePairs().begin();
		for(; it != aModel.getBasePairs().end(); ++ it)
		{
			bool bAdjacent = false;
			bAdjacent = areAdjacent(
				(*it).fResId, 
				(*it).rResId, 
				associationMap);
				
			if(!bAdjacent)
			{
				mPairs.push_back(*it);
			}			
		}
	}
	
	bool AnnotationTertiaryPairs::areAdjacent(
		const mccore::ResId& aResId1, 
		const mccore::ResId& aResId2, 
		const struct_association_map& associationMap) const
	{
		bool bAdjacent = false;
		
		std::pair<
			struct_association_map::const_iterator, 
			struct_association_map::const_iterator> range1;
			
		std::pair<
			struct_association_map::const_iterator, 
			struct_association_map::const_iterator> range2;
			
		range1 = associationMap.equal_range(aResId1);
		range2 = associationMap.equal_range(aResId2);
		
		struct_association_map::const_iterator it1;
		struct_association_map::const_iterator it2;
		
		for(it1 = range1.first; it1 != range1.second && !bAdjacent; ++ it1)
		{
			const SecondaryStructure* pStruct1 = it1->second;
			for(it2 = range2.first; it2 != range2.second && !bAdjacent; ++ it2)
			{
				bAdjacent |= pStruct1->isAdjacent(*it2->second);
			}		
		}
		
		return bAdjacent;
	}
	
	void AnnotationTertiaryPairs::addStemsAssociations(
		const AnnotateModel& aModel, 
		struct_association_map& associationMap	)
	{
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>(AnnotationStems().provides());
		
		if(NULL != pAnnotStems)
		{
			GraphModel::const_iterator it = aModel.begin();
			for(;it != aModel.end(); ++ it)
			{
				const ResId resId = (*it).getResId();
				
				// Add all the stem associations
				std::vector<Stem>::const_iterator stemIt;
				for(stemIt = pAnnotStems->getStems().begin(); 
					stemIt != pAnnotStems->getStems().end(); 
					++stemIt)
				{
					if(stemIt->contains(resId))
					{
						struct_association	mapping(resId, &(*stemIt));
						associationMap.insert(mapping);
					}
				}			
			}
		}
	}
	
	void AnnotationTertiaryPairs::addLoopsAssociations(
		const AnnotateModel& aModel, 
		struct_association_map& associationMap	)
	{
		// Add the loop associations
		const AnnotationLoops* pAnnotLoops = NULL;
		pAnnotLoops = aModel.getAnnotation<AnnotationLoops>("Loops");

		if(NULL != pAnnotLoops)
		{		
			GraphModel::const_iterator it = aModel.begin();
			for(;it != aModel.end(); ++ it)
			{
				const ResId resId = (*it).getResId();
				
				if(NULL != pAnnotLoops)
				{
					std::vector<Loop>::const_iterator loopIt;
					loopIt = pAnnotLoops->getLoops().begin();
					for(; loopIt != pAnnotLoops->getLoops().end(); ++ loopIt)
					{
						if(loopIt->contains(resId))
						{
							struct_association mapping(resId, &(*loopIt));
							associationMap.insert(mapping);
						}
					}
				}			
			}
		}
	}
	
	std::string AnnotationTertiaryPairs::output() const
	{
		std::ostringstream oss;
		std::vector<BasePair>::const_iterator it;
		for(it = mPairs.begin(); it != mPairs.end(); ++ it)
		{
			oss << it->fResId << "-" << it->rResId << std::endl;
		}
		return oss.str();
	}
}
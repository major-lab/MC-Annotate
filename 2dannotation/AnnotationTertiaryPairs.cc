#include "AnnotationTertiaryPairs.h"
#include "AnnotateModel.h"
#include "AnnotationInteractions.h"
#include "AnnotationStems.h"
#include "AnnotationLoops.h"
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationTertiaryPairs::mstrAnnotationName = "Tertiary Pairs";
	
	// Methods
	AnnotationTertiaryPairs::AnnotationTertiaryPairs()
	{
		// Requires
		addRequirement<AnnotationInteractions>();
		addRequirement<AnnotationStems>();
		addRequirement<AnnotationLoops>();
	}
	
	AnnotationTertiaryPairs::~AnnotationTertiaryPairs()
	{
		clear();
	}
	
	void AnnotationTertiaryPairs::clear()
	{
		mPairs.clear();
	}
	
	const std::set< BasePair >& AnnotationTertiaryPairs::getPairs() const
	{
		return mPairs;
	}
	
	void AnnotationTertiaryPairs::update(AnnotateModel& aModel)
	{
		struct_association_map associationMap;
		
		const AnnotationInteractions* pAInteractions = NULL;
		pAInteractions = aModel.getAnnotation<AnnotationInteractions>();
		
		if(NULL != pAInteractions)
		{
			// Add the stem associations
			addStemsAssociations(aModel, associationMap);
			
			// Add the loop associations
			addLoopsAssociations(aModel, associationMap);
			
			std::vector<BasePair>::const_iterator it = pAInteractions->getPairs().begin();
			for(; it != pAInteractions->getPairs().end(); ++ it)
			{
				// Filter out pairs between different chains
				if(it->fResId.getChainId() == it->rResId.getChainId())
				{
					// Filter out adjacent pairs
					bool bAdjacent = false;
					bAdjacent = areAdjacent(it->fResId, it->rResId, associationMap);
					if(!bAdjacent)
					{
						// Filter out pairs inside a loop
						bool bSameLoop = false;
						bSameLoop = areSameLoop(aModel, it->fResId, it->rResId);
						if(!bSameLoop)
						{
							mPairs.insert(*it);
						}
					}
				}
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
				bAdjacent = pStruct1->isAdjacent(*it2->second);
			}		
		}
		
		return bAdjacent;
	}
	
	bool AnnotationTertiaryPairs::areSameLoop(
		const AnnotateModel& aModel,
		const mccore::ResId& aResId1, 
		const mccore::ResId& aResId2) const
	{
		bool bSameLoop = false;
		
		const AnnotationLoops* pAnnotLoops = 
			aModel.getAnnotation<AnnotationLoops>();
			
		if(NULL != pAnnotLoops)
		{
			std::vector< Loop >::const_iterator it;
			for(it = pAnnotLoops->getLoops().begin(); 
				it != pAnnotLoops->getLoops().end() && !bSameLoop; 
				++ it)
			{
				std::set<mccore::ResId> resids = it->getResIds();
				if(	resids.end() != resids.find(aResId1) 
					&& resids.end() != resids.find(aResId2))
				{
					bSameLoop = true;
				}
			}
		}		
		return bSameLoop;
	}
	
	void AnnotationTertiaryPairs::addStemsAssociations(
		const AnnotateModel& aModel, 
		struct_association_map& associationMap	)
	{
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>();
		
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
		pAnnotLoops = aModel.getAnnotation<AnnotationLoops>();

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
		std::set<BasePair>::const_iterator it;
		for(it = mPairs.begin(); it != mPairs.end(); ++ it)
		{
			oss << it->fResId << "-" << it->rResId << std::endl;
		}
		return oss.str();
	}
}
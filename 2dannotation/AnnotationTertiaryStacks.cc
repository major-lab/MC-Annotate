#include "AnnotationTertiaryStacks.h"
#include "AnnotateModel.h"
#include "AnnotationInteractions.h"
#include "AnnotationStems.h"
#include "AnnotationLoops.h"
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationTertiaryStacks::mstrAnnotationName = "Tertiary Stacks";
	
	// Methods
	AnnotationTertiaryStacks::AnnotationTertiaryStacks()
	{
		// Requires
		addRequirement<AnnotationInteractions>();
		addRequirement<AnnotationStems>();
		addRequirement<AnnotationLoops>();
	}
	
	AnnotationTertiaryStacks::~AnnotationTertiaryStacks()
	{
		clear();
	}
	
	void AnnotationTertiaryStacks::clear()
	{
		mStacks.clear();
	}
	
	const std::set< BaseStack >& AnnotationTertiaryStacks::getStacks() const
	{
		return mStacks;
	}
	
	// TODO : Secondary structure association should be outside of these to be 
	// computed only once
	void AnnotationTertiaryStacks::update(const AnnotateModel& aModel)
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
			
			std::vector<BaseStack>::const_iterator it = pAInteractions->getStacks().begin();
			for(; it != pAInteractions->getStacks().end(); ++ it)
			{
				bool bAdjacent = false;
				bAdjacent = areAdjacent(
					(*it).fResId, 
					(*it).rResId, 
					associationMap);
					
				if(!bAdjacent)
				{
					mStacks.insert(*it);
				}			
			}
		}
	}
	
	bool AnnotationTertiaryStacks::areAdjacent(
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
	
	void AnnotationTertiaryStacks::addStemsAssociations(
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
	
	void AnnotationTertiaryStacks::addLoopsAssociations(
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
	
	std::string AnnotationTertiaryStacks::output() const
	{
		std::ostringstream oss;
		std::set<BaseStack>::const_iterator it;
		for(it = mStacks.begin(); it != mStacks.end(); ++ it)
		{
			oss << it->fResId << "-" << it->rResId << std::endl;
		}
		return oss.str();
	}
}
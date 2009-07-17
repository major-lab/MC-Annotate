#include "AnnotationResSecondaryStructures.h"
#include "AnnotateModel.h"
#include "AnnotationLoops.h"
#include "AnnotationStems.h"
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationResSecondaryStructures::mstrAnnotationName = "ResSecondaryStructures";
	
	// Methods	
	AnnotationResSecondaryStructures::AnnotationResSecondaryStructures() 
	{
		addRequirement<AnnotationLoops>();
		addRequirement<AnnotationStems>();
	}
	
	AnnotationResSecondaryStructures::~AnnotationResSecondaryStructures() 
	{
		clear();
	}
	
	void AnnotationResSecondaryStructures::clear()
	{
		mMapping.clear();
	}
		
	void AnnotationResSecondaryStructures::update(const AnnotateModel& aModel)
	{
		// This model contains no stems, there is only one open loop
		const AnnotationLoops* pAnnotLoops = aModel.getAnnotation<AnnotationLoops>();
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>();
		
		if(NULL != pAnnotLoops && NULL != pAnnotStems)
		{
			const SecondaryStructure* pStruct = NULL;
			mccore::GraphModel::const_iterator itRes = aModel.begin();
			for(; itRes != aModel.end(); ++ itRes)
			{
				pStruct = NULL;
				mccore::ResId resId = itRes->getResId();
				std::vector< Stem >::const_iterator itStem;
				for(itStem = pAnnotStems->getStems().begin(); 
					itStem != pAnnotStems->getStems().end() && NULL == pStruct; 
					++itStem)
				{
					if(itStem->contains(resId))
					{
						pStruct = &(*itStem);
					}
				}
				
				std::vector< Loop >::const_iterator itLoop;
				for(itLoop = pAnnotLoops->getLoops().begin(); 
					itLoop != pAnnotLoops->getLoops().end() && NULL == pStruct; 
					++itLoop)
				{
					if(itLoop->contains(resId))
					{
						pStruct = &(*itLoop);
					}
				}
				
				if(NULL == pStruct)
				{
					mccore::gOut(0) << "Residue " << resId;
					mccore::gOut(0) << " is associated with no structure";
					mccore::gOut(0) << std::endl;
					
				}
				else
				{
					mMapping[resId] = pStruct;
				}
			}			
		}
	}
	
	std::string AnnotationResSecondaryStructures::output() const
	{
		std::ostringstream oss;
		return oss.str();
	}
}
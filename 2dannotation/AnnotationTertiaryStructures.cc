#include "AnnotationTertiaryStructures.h"

#include <sstream>

#include "mccore/Pdbstream.h"

#include "AnnotateModel.h"
#include "AnnotationStems.h"
#include "AnnotationLoops.h"
#include "AnnotationResSecondaryStructures.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotationTertiaryStacks.h"
#include "AnnotationTertiaryCycles.h"

namespace annotate
{
	std::string AnnotationTertiaryStructures::mstrAnnotationName = "Tertiary Structures";
	
	AnnotationTertiaryStructures::AnnotationTertiaryStructures()
	{
		addRequirement<AnnotationStems>();
		addRequirement<AnnotationLoops>();
		addRequirement<AnnotationTertiaryPairs>();
		addRequirement<AnnotationTertiaryStacks>();
		addRequirement<AnnotationTertiaryCycles>();
		addRequirement<AnnotationResSecondaryStructures>();
	}
	
	AnnotationTertiaryStructures::~AnnotationTertiaryStructures()
	{
		clear();
	}
	
	void AnnotationTertiaryStructures::clear()
	{
		mSinglePairs.clear();
		mSingleStacks.clear();
		mStructures.clear();
	}
	
	void AnnotationTertiaryStructures::update(AnnotateModel& aModel)
	{
		clear();
				
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>();
		
		// Extract the single pairs
		const AnnotationTertiaryCycles* pAnnotCycles = 
			aModel.getAnnotation<AnnotationTertiaryCycles>();
		const AnnotationTertiaryPairs* pAnnotPairs = 
			aModel.getAnnotation<AnnotationTertiaryPairs>();
		const AnnotationTertiaryStacks* pAnnotStacks = 
			aModel.getAnnotation<AnnotationTertiaryStacks>();
			
		if(	NULL != pAnnotStems && NULL != pAnnotCycles 
			&& NULL != pAnnotPairs && NULL != pAnnotStacks)
		{
			std::list<Cycle> cycles;
			std::set<Cycle>::const_iterator itSet = pAnnotCycles->getCycles().begin();
			for(; itSet != pAnnotCycles->getCycles().end(); ++itSet)
			{
				cycles.push_back(*itSet);
			}			
			
			// Get all the interactions from the cycles
			std::set<BaseInteraction> cycleInteractions;
			std::list< Cycle >::iterator itCycle;
			for(itCycle = cycles.begin(); itCycle != cycles.end(); ++itCycle)
			{
				pAnnotCycles->getPairs(*itCycle, cycleInteractions);			
			}
			
			// Add the intersection of the tertiary pair and the cycle's interaction 
			// to the set of single pairs structure
			mSinglePairs = difference(pAnnotPairs->getPairs(), cycleInteractions);
			
			// Add the intersection of the tertiary stacks and the cycle's interaction 
			// to the set of single pairs structure
			mSingleStacks = difference(pAnnotStacks->getStacks(), cycleInteractions);
			
			while(!cycles.empty())
			{
				TertiaryStructure structure;
				structure.addCycle(cycles.front());
				cycles.pop_front();
				itCycle = cycles.begin();
				while(itCycle != cycles.end())
				{
					if(itCycle->shareInteractions(structure.getCycles().front()))
					{
						structure.addCycle(*itCycle);
						itCycle = cycles.erase(itCycle);
					}
					else
					{
						itCycle ++;
					}
				}
				mStructures.push_back(structure);
			}
		}
		
		// Compute models from all the structures
		computeModels(aModel);
	}
	
	std::set<BasePair> AnnotationTertiaryStructures::difference(
		const std::set<BasePair>& a3DPair, 
		const std::set<BaseInteraction>& aCyclePairs) const
	{
		std::set<BasePair> pairs;
		std::set<BasePair>::const_iterator first1 = a3DPair.begin();
		std::set<BasePair>::const_iterator last1 = a3DPair.end();
		std::set<BaseInteraction>::const_iterator first2 = aCyclePairs.begin();
		std::set<BaseInteraction>::const_iterator last2 = aCyclePairs.end();		
		while (first1!=last1 && first2!=last2)
  		{
  			const BaseInteraction* pLeft = &(*first1);
  			const BaseInteraction* pRight = &(*first2);
  			
			if (*pLeft<*pRight)
			{
				pairs.insert(*first1);
				 ++first1;
			}
			else if (*pRight<*pLeft)
			{
				 ++first2;
			}
			else
			{
				++first1;
				++first2;
			}
		}				
		return pairs;		
	}
	
	std::set<BaseStack> AnnotationTertiaryStructures::difference(
		const std::set<BaseStack>& a3DStack, 
		const std::set<BaseInteraction>& aCyclePairs) const
	{
		std::set<BaseStack> pairs;
		std::set<BaseStack>::const_iterator first1 = a3DStack.begin();
		std::set<BaseStack>::const_iterator last1 = a3DStack.end();
		std::set<BaseInteraction>::const_iterator first2 = aCyclePairs.begin();
		std::set<BaseInteraction>::const_iterator last2 = aCyclePairs.end();		
		while (first1!=last1 && first2!=last2)
  		{
  			const BaseInteraction* pLeft = &(*first1);
  			const BaseInteraction* pRight = &(*first2);
  			
			if (*pLeft<*pRight)
			{
				pairs.insert(*first1);
				 ++first1;
			}
			else if (*pRight<*pLeft)
			{
				 ++first2;
			}
			else
			{
				++first1;
				++first2;
			}
		}				
		return pairs;		
	}
	
	std::string AnnotationTertiaryStructures::output() const
	{
		std::ostringstream oss;
		
		oss << "Single pair tertiary structure ----------------------------------" << std::endl;
		std::set<BasePair>::const_iterator itPair;
		for(itPair = mSinglePairs.begin(); itPair != mSinglePairs.end(); ++ itPair)
		{
			oss << itPair->fResId << "-" << itPair->rResId << std::endl;
		}
		
		oss << "Single stack tertiary structure ---------------------------------" << std::endl;
		std::set<BaseStack>::const_iterator itStack;
		for(itStack = mSingleStacks.begin(); itStack != mSingleStacks.end(); ++ itStack)
		{
			oss << itStack->fResId << "-" << itStack->rResId << std::endl;
		}
		
		oss << "Multiple cycle tertiary structure -------------------------------" << std::endl;
		std::list<TertiaryStructure>::const_iterator itStruct;
		for(itStruct = mStructures.begin(); itStruct != mStructures.end(); ++itStruct)
		{
			oss << itStruct->name();
			oss << std::endl;
		}
		
		return oss.str();
	}
	
	std::string AnnotationTertiaryStructures::structureName(
		const AnnotateModel& aModel, 
		const TertiaryStructure& aStructure) const
	{
		std::ostringstream oss;
		
		const AnnotationResSecondaryStructures* pAnnotSecondaryStruct = 
			aModel.getAnnotation<AnnotationResSecondaryStructures>();
		
		// Add the name of the model where it comes from
		oss << aModel.name();
		
		// Add the residues id in the model
		std::set<mccore::ResId>::const_iterator it;
		for(it = aStructure.getResidues().begin();
			it != aStructure.getResidues().end(); 
			++ it)
		{
			if(it == aStructure.getResidues().begin())
			{
				oss << "-";
			}
			else
			{
				oss << ".";
			}
			oss << *it;
		}
		
		// Adds the structure
		std::set<mccore::ResId>::const_iterator itPrev = aStructure.getResidues().end();
		for(it = aStructure.getResidues().begin();
			it != aStructure.getResidues().end(); 
			++ it)
		{
			if(it == aStructure.getResidues().begin())
			{
				oss << "-";
			}
			else if(!itPrev->areContiguous(*it))
			{
				oss << "_";
			}
			oss << Pdbstream::stringifyResidueType (aModel.find(*it)->getType());
			itPrev = it;
		}
		
		// Add the associated secondary structure
		itPrev = aStructure.getResidues().end();
		const SecondaryStructure* pPrevStruct = NULL;
		for(it = aStructure.getResidues().begin();
			it != aStructure.getResidues().end(); 
			++ it)
		{
			if(it == aStructure.getResidues().begin())
			{
				oss << "-";
			}
			else if(!itPrev->areContiguous(*it))
			{
				oss << "_";
				pPrevStruct = NULL;
			}
			
			
			if(NULL != pAnnotSecondaryStruct)
			{
				std::map<mccore::ResId, const SecondaryStructure*>::const_iterator itMap;
				itMap = pAnnotSecondaryStruct->getMapping().find(*it);
				if(pAnnotSecondaryStruct->getMapping().end() != itMap)
				{
					const SecondaryStructure* pStruct = itMap->second;
					if(pPrevStruct != pStruct)
					{
						if(pPrevStruct != NULL)
						{
							oss << ".";
						}		
						const Stem* pStem = dynamic_cast<const Stem*>(pStruct);
						if(NULL != pStem)
						{
							oss << "stem";
						}
						else
						{
							const Loop* pLoop = dynamic_cast<const Loop*>(pStruct);
							if(NULL != pLoop)
							{
								oss << pLoop->describe() << "loop";
							}
							else
							{
								oss << "unknown";
							}
						}
						pPrevStruct = pStruct;
					}
				}
			}
			itPrev = it;
		}
		
		return oss.str();
	}
	
	void AnnotationTertiaryStructures::computeModels(const AnnotateModel& aModel)
	{
		std::list<TertiaryStructure>::iterator itStruct;
		for(itStruct = mStructures.begin(); itStruct != mStructures.end(); ++itStruct)
		{
			itStruct->update(aModel);
			itStruct->name(structureName(aModel, *itStruct));			
		}
	}
}
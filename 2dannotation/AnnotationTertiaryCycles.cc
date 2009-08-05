#include "AnnotationTertiaryCycles.h"
#include "AnnotationCycles.h"
#include "AnnotationInteractions.h"
#include "AnnotationLoops.h"
#include "AnnotationResSecondaryStructures.h"
#include "AnnotationTertiaryPairs.h"
#include "AnnotationTertiaryStacks.h"
#include "AnnotateModel.h"
#include "AlgorithmExtra.h"
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationTertiaryCycles::mstrAnnotationName = "Tertiary Cycles";
	
	// Methods	
	AnnotationTertiaryCycles::AnnotationTertiaryCycles(unsigned int auiMaxCycleSize)
	{
		addRequirement<AnnotationCycles>();
		addRequirement<AnnotationLoops>();
		addRequirement<AnnotationTertiaryPairs>();
		addRequirement<AnnotationTertiaryStacks>();
		addRequirement<AnnotationResSecondaryStructures>();
		muiMaxCycleSize = auiMaxCycleSize;
	}
	
	AnnotationTertiaryCycles::~AnnotationTertiaryCycles()
	{
		clear();
	}
	
	void AnnotationTertiaryCycles::clear()
	{
		mCycles.clear();
		mConnects.clear();
	}
		
	void AnnotationTertiaryCycles::update(const AnnotateModel& aModel)
	{
		clear();
		
		const AnnotationCycles* pAnnCycles = 
			aModel.getAnnotation<AnnotationCycles>();
		
		if(	NULL != pAnnCycles)
		{
			// Copy all the cycles from the annotation of cycles
			std::list<Cycle> cycles;
			cycles.insert(
				cycles.begin(), 
				pAnnCycles->getCycles().begin(), 
				pAnnCycles->getCycles().end());
				
			// Create cycles for the loops
			updateCycleLoops(aModel);
			
			// From the list of potential tertiary cycles, 
			// keep only those with tertiary interactions
			filterOutAdjacentStructures(aModel, cycles);
			
			// From the list of potential tertiary cycles, 
			// remove those that are subparts of a loop
			filterOutLoops(cycles);
			
			// Keep the remaining loops as tertiary
			mCycles.insert(cycles.begin(), cycles.end());
		}
		updateConnections(aModel);
	}
	
	void AnnotationTertiaryCycles::filterOutAdjacentStructures(
		const AnnotateModel& aModel, 
		std::list<Cycle>& aCycles) const
	{
		const AnnotationTertiaryPairs* pAnnTPairs = NULL;
		const AnnotationTertiaryStacks* pAnnTStacks = NULL;
		pAnnTPairs = aModel.getAnnotation<AnnotationTertiaryPairs>();
		pAnnTStacks = aModel.getAnnotation<AnnotationTertiaryStacks>();
		
		if(NULL != pAnnTPairs && NULL != pAnnTStacks)
		{
			// From all the cycles, keep only those that are non-adjacent
			std::set<BaseInteraction> cyclePairs;
			std::list<Cycle>::iterator it = aCycles.begin();
			while(it != aCycles.end())
			{
				bool bTertiary = false;
				unsigned int uiSize = it->getResidues().size();
				if((0 == muiMaxCycleSize || uiSize < muiMaxCycleSize) 
					&& it->isSingleChain())
				{
					// Check if the cycle is tertiary
					getPairs(*it, cyclePairs);
				
					bTertiary = isTertiary(
						cyclePairs, 
						pAnnTPairs->getPairs());
					if(!bTertiary)
					{
						bTertiary = isTertiary(
							cyclePairs, 
							pAnnTStacks->getStacks());
					}
					cyclePairs.clear();
				}
				
				if(bTertiary)
				{
					++it;
				}
				else
				{
					it = aCycles.erase(it);
				}
			}
		}
	}
	
	void AnnotationTertiaryCycles::filterOutLoops(
		std::list<Cycle>& aCycles) const
	{		
		std::list<Cycle>::iterator it = aCycles.begin();
					
		while(it != aCycles.end())
		{
			bool bLoopPart = false;
			std::set<Cycle>::const_iterator itLoop;
			for(itLoop = mLoops.begin(); 
				itLoop != mLoops.end() && !bLoopPart; 
				++itLoop)
			{
				std::set<BaseInteraction> interactions = 
					it->getBaseInteractions();
				unsigned int uiDiffSize = SetDifference(
					interactions, 
					itLoop->getBaseInteractions()).size();
				bLoopPart = (0 == uiDiffSize);
			}
			if(bLoopPart)
			{
				it = aCycles.erase(it);
			}
			else
			{
				++ it;
			}		
		}
	}
	
	void AnnotationTertiaryCycles::updateCycleLoops(const AnnotateModel& aModel)
	{
		mLoops.clear();
		
		const AnnotationLoops* pAnnLoops = 
			aModel.getAnnotation<AnnotationLoops>();
		
		if(NULL != pAnnLoops)
		{
			unsigned char ucRelMask = aModel.relationMask();
			std::vector<Loop> loops = pAnnLoops->getLoops();
			std::vector<Loop>::const_iterator itLoop;
			for(itLoop = loops.begin(); itLoop != loops.end(); ++ itLoop)
			{
				std::set<BaseInteraction> interactions = itLoop->getBaseInteractions();
				// A cycle requires at least 3 interactions
				if(2 < interactions.size()) 
				{
					Cycle cycle(aModel, interactions, ucRelMask);
					mLoops.insert(cycle);
				}
			}
		}
	}
	
	void AnnotationTertiaryCycles::updateConnections(const AnnotateModel& aModel)
	{
		const AnnotationCycles* pACycles = 
			aModel.getAnnotation<AnnotationCycles>();
		
		if(NULL != pACycles)
		{
			std::list<Cycle> cycles;
			
				
			std::set<Cycle> adjacents;
			adjacents = SetDifference<Cycle>(pACycles->getCycles(), mCycles);
			cycles.insert(cycles.begin(), adjacents.begin(), adjacents.end());
			filterOutLoops(cycles);
			std::set<Cycle>::const_iterator it;
			for(it = mCycles.begin(); it != mCycles.end(); ++ it)
			{
				std::set<Cycle> connections;
				
				// Add the 2D cycles associated
				std::list<Cycle>::const_iterator adjIt;
				for(adjIt = cycles.begin(); adjIt != cycles.end(); ++adjIt)
				{
					if(it->shareInteractions(*adjIt))
					{
						connections.insert(*adjIt);
					}
				}
				
				// Add the loops cycles
				std::set<Cycle>::const_iterator itCycle;
				for(itCycle = mLoops.begin(); itCycle != mLoops.end(); ++itCycle)
				{
					if(it->shareInteractions(*itCycle))
					{
						connections.insert(*itCycle);
					}
				}
								
				mConnects.insert(std::pair<Cycle, std::set<Cycle> >(*it, connections));
				std::list< std::list<mccore::ResId> > linkConnects = updateLinkerConnections(aModel, *it);
				mConnectsLinkers.insert(std::pair<Cycle, linkers_connect>(*it, linkConnects));
			}
		}
	}
	
	std::list< std::list<mccore::ResId> > 
	AnnotationTertiaryCycles::updateLinkerConnections(
		const AnnotateModel& aModel, 
		const Cycle& aCycle) const
	{
		std::list< std::list<mccore::ResId> > connections;
		
		const AnnotationResSecondaryStructures* pAResIdStruct = 
			aModel.getAnnotation<AnnotationResSecondaryStructures>();
			
		if(NULL != pAResIdStruct)
		{
			const Loop* pPrevLoop = NULL;
			std::list<mccore::ResId> strand;
			std::list<mccore::ResId>::const_iterator it;
			for(it = aCycle.getResidues().begin(); 
				it != aCycle.getResidues().end(); 
				++it)
			{
				std::map<mccore::ResId, const SecondaryStructure*>::const_iterator mapIt;
				mapIt = pAResIdStruct->getMapping().find(*it);
				if(mapIt != pAResIdStruct->getMapping().end())
				{
					// Is this an open loop ?
					const Loop* pLoop = dynamic_cast<const Loop*>(mapIt->second);
					
					if(NULL != pLoop && 0 == pLoop->describe().compare("open"))
					{
						// This is an open loop
						if(0 != strand.size() && (pLoop != pPrevLoop))
						{
							connections.push_back(strand);
							strand.clear();							
						}
						strand.push_back(*it);
						pPrevLoop = pLoop;
					}
				}
			}
			if(0 < strand.size())
			{
				connections.push_back(strand);
				strand.clear();
			}
		}
		return connections;
	}
	
	std::string AnnotationTertiaryCycles::output() const
	{
		std::ostringstream oss;
		std::set<Cycle>::const_iterator itCycle;
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
	
	const std::set<Cycle>& AnnotationTertiaryCycles::getConnections(const Cycle& aCycle) const
	{
		std::map< Cycle, std::set <Cycle> >::const_iterator it = mConnects.find(aCycle);
		if(it == mConnects.end())
		{
			throw mccore::NoSuchElementException("Cycle not found", __FILE__, __LINE__);
		}
		return it->second;
	}
	
	const std::list<std::list<mccore::ResId> >& AnnotationTertiaryCycles::getLinkerConnections(
		const Cycle& aCycle) const
	{
		std::map< Cycle, linkers_connect >::const_iterator it = mConnectsLinkers.find(aCycle);
		if(it == mConnectsLinkers.end())
		{
			throw mccore::NoSuchElementException("Linkers not found", __FILE__, __LINE__);
		}
		return it->second;
	}
}
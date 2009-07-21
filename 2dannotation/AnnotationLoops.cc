#include "AnnotationLoops.h"
#include "AnnotateModel.h"
#include "AnnotationLinkers.h"
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationLoops::mstrAnnotationName = "Loops";
	
	// Methods	
	AnnotationLoops::AnnotationLoops() 
	{
		addRequirement<AnnotationLinkers>();
	}
	
	AnnotationLoops::~AnnotationLoops() 
	{
		clear();
	}
	
	void AnnotationLoops::clear()
	{
		mLoops.clear();
	}
	
	const std::vector< Loop >& AnnotationLoops::getLoops() const
	{
		return mLoops;
	}
	
	void AnnotationLoops::update(const AnnotateModel& aModel)
	{
		// This model contains no stems, there is only one open loop
		const AnnotationLinkers* pAnnotLinkers = NULL;
		pAnnotLinkers = aModel.getAnnotation<AnnotationLinkers>();
		
		if(NULL != pAnnotLinkers)
		{
			/*
			std::list<Loop> potentialLoops;
			std::vector<Linker>::const_iterator itLinker;
			mccore::gOut(0) << pAnnotLinkers->getLinkers().size() << " linkers" << std::endl;
			for(itLinker = pAnnotLinkers->getLinkers().begin();
				itLinker != pAnnotLinkers->getLinkers().end(); 
				++ itLinker)
			{
				potentialLoops.push_back(Loop(*itLinker));
			}
			mccore::gOut(0) << "Identified " << potentialLoops.size() << " potential loops" << std::endl;
			
			while(!potentialLoops.empty())
			{
				if(potentialLoops.front().complete())
				{
					mLoops.push_back(potentialLoops.front());
					potentialLoops.pop_front();
					mccore::gOut(0) << "Loop identified" << std::endl;
					mccore::gOut(0) << potentialLoops.size() << " potential loops left" << std::endl;
					
				}else
				{
					Loop workLoop = potentialLoops.front();
					potentialLoops.pop_front();

					// Find other loop to connect this one to
					std::list<Loop>::iterator it = potentialLoops.begin();
					while(it != potentialLoops.end())
					{
						if(workLoop.connects(*it))
						{
							workLoop.append(*it);
							it = potentialLoops.erase(it);
						}
						else
						{
							it ++;
						}
					}
					potentialLoops.push_back(workLoop);
					mccore::gOut(0) << potentialLoops.size() << " potential loops left" << std::endl;
				}
			}*/

			std::map<ResId, const Linker*> residueLinkerMap;
			residueLinkerMap = getResidueLinkerMap(*pAnnotLinkers);
			if(residueLinkerMap.empty() && 0 < aModel.size())
			{
				// TODO : Assert that there is only one linker
				if(0 < pAnnotLinkers->getLinkers().size())
				{
					// Add look if there is one
					std::vector<Linker> potentialLoop;
					potentialLoop.push_back(pAnnotLinkers->getLinkers().front());
					mLoops.push_back(Loop(potentialLoop));
				}
			}
			else
			{
				// Find all the loops contained in the linker map
				findLoops(residueLinkerMap);
			}
		}
	}
	
	void AnnotationLoops::findLoops(
		std::map<ResId, const Linker*> residueLinkerMap)
	{
		std::vector< Loop > openLoops;
		while(!residueLinkerMap.empty())
		{
			const Linker* pFirstLinker = NULL;
			std::vector<Linker> potentialLoop;
			
			// Pick the first linker
			std::pair<ResId, const Linker*> mapping = *(residueLinkerMap.begin());
			
			// Navigate until next is stem is NULL or we loop
			pFirstLinker = mapping.second;
			potentialLoop.push_back(*pFirstLinker);
			
			// Remove the linker from the map
			removeLinker(residueLinkerMap, *pFirstLinker);
			
			Linker currentLinker = nextLinker(
				*pFirstLinker, 
				residueLinkerMap);
			
			while(!currentLinker.isEmpty() && currentLinker != *pFirstLinker)
			{
				potentialLoop.push_back(currentLinker);
				
				// Remove the linker from the map
				removeLinker(residueLinkerMap, currentLinker);
				
				// Get the next strands
				currentLinker = nextLinker(
					currentLinker, 
					residueLinkerMap);				
			}
			
			if(0 < potentialLoop.size())
			{
				const StemConnection* pStart = &potentialLoop.front().getStart();
				const StemConnection* pEnd = &potentialLoop.back().getEnd();
				if(	pStart->isValid() 
					&& pEnd->isValid() 
					&& (&pStart->getStem() == &pEnd->getStem()))
				{
					mLoops.push_back(Loop(potentialLoop));
				}
				else
				{
					openLoops.push_back(Loop(potentialLoop));
				}
			}
		}
		// Find the open loops from the remaining linkers
		findOpenLoops(openLoops);
	}
	
	void AnnotationLoops::findOpenLoops(std::vector< Loop >& aOpenLoops)
	{
		while(!aOpenLoops.empty())
		{
			std::vector<Loop>::iterator it = aOpenLoops.begin();
			Loop openLoop = *it;
			aOpenLoops.erase(it);
			const Linker* pFrontLinker = &openLoop.getLinkers().front();
			const Linker* pBackLinker = &openLoop.getLinkers().back();
			
			if(pFrontLinker->getStart().isValid())
			{
				// Starts by a stem, check if another open loop precedes it
				const ResId resId = pFrontLinker->getStart().nextId();
				std::vector<Loop>::iterator potentialIt = 
					getLoopStartingBy(resId, aOpenLoops);
				
				if(potentialIt != aOpenLoops.end())
				{
					(*potentialIt).reverse();
					(*potentialIt).append(openLoop);
				}else
				{
					potentialIt = getLoopEndingBy(resId, aOpenLoops);
					if(potentialIt != aOpenLoops.end())
					{
						(*potentialIt).append(openLoop);						
					}
					else
					{
						mLoops.push_back(openLoop);
					}
				}				
			}else if(pBackLinker->getEnd().isValid())
			{
				// Ends by a stem, check if another open loop continues it
				const ResId resId = pBackLinker->getEnd().nextId();
				std::vector<Loop>::iterator potentialIt = 
					getLoopStartingBy(resId, aOpenLoops);
				
				if(potentialIt != aOpenLoops.end())
				{
					(*potentialIt).append(openLoop);
				}else
				{
					potentialIt = getLoopEndingBy(resId, aOpenLoops);
					if(potentialIt != aOpenLoops.end())
					{
						(*potentialIt).reverse();
						(*potentialIt).append(openLoop);						
					}
					else
					{
						mLoops.push_back(openLoop);
					}
				}
			}
			else
			{
				// The open loop is complete, add it to the loops
				mLoops.push_back(openLoop);
			}
		}		
	}
	
	std::vector<Loop>::iterator 
	AnnotationLoops::getLoopStartingBy(
		const ResId& aResId, 
		std::vector<Loop>& aLoops) const
	{
		std::vector<Loop>::iterator it = aLoops.begin();
		for(; it != aLoops.end(); ++ it)
		{
			if((*it).getLinkers().front().getStart().isValid() 
				&& aResId == (*it).getLinkers().front().getStart().getResidue())
			{
				break;				
			}
		}
		return it;		
	}
	
	std::vector<Loop>::iterator 
	AnnotationLoops::getLoopEndingBy(
		const ResId& aResId, 
		std::vector<Loop>& aLoops) const
	{
		std::vector<Loop>::iterator it = aLoops.begin();
		for(; it != aLoops.end(); ++ it)
		{
			if((*it).getLinkers().back().getEnd().isValid()
				&& aResId == (*it).getLinkers().back().getEnd().getResidue())
			{
				break;				
			}
		}
		return it;		
	}
	
	std::map<mccore::ResId, const Linker*>
	AnnotationLoops::getResidueLinkerMap(const AnnotationLinkers& aAnnotLinkers) const
	{
		std::map<mccore::ResId, const Linker*> residueLinkerMap;
		std::vector<Linker>::const_iterator it;
		for(it = aAnnotLinkers.getLinkers().begin();
			it != aAnnotLinkers.getLinkers().end(); 
			++it)
		{
			if((*it).getStart().isValid())
			{
				const StemConnection* pConnect = &(*it).getStart();
				ResId resId = pConnect->getResidue();
				std::pair<ResId, const Linker*> mapping(resId, &(*it));
				residueLinkerMap.insert(mapping);	
			}
			
			if((*it).getEnd().isValid())
			{
				const StemConnection* pConnection = &(*it).getEnd();
				ResId resId = pConnection->getResidue();
				std::pair<ResId, const Linker*> mapping(resId, &(*it));
				residueLinkerMap.insert(mapping);
			}
		}
		return residueLinkerMap;
	}
	
	void AnnotationLoops::removeLinker(
		std::map<mccore::ResId, const Linker*>& aResidueLinkerMap,
		const Linker& aLinker) const
	{
		if(aLinker.getEnd().isValid())
		{
			ResId resId = aLinker.getEnd().getResidue();
			aResidueLinkerMap.erase(resId);
		}
		
		if(aLinker.getStart().isValid())
		{
			ResId resId = aLinker.getStart().getResidue();
			aResidueLinkerMap.erase(resId);
		}
	}
	
	Linker AnnotationLoops::nextLinker(
		const Linker& aLinker,
		const std::map<mccore::ResId, const Linker*>& aResidueLinkerMap) const
	{
		Linker nextLinker;
		const Stem* pNextStem = &aLinker.getEnd().getStem();
		if(NULL != pNextStem)
		{
			ResId resId = aLinker.getEnd().nextId();
			std::map<mccore::ResId, const Linker*>::const_iterator it;
			it = aResidueLinkerMap.find(resId);
			if(it != aResidueLinkerMap.end())
			{
				const Linker* pLinker = (*it).second;
				nextLinker = *pLinker;
				
				if(nextLinker.getEnd().isValid())
				{
					ResId endId = nextLinker.getEnd().getResidue();
					if(it->first == endId)
					{
						nextLinker.reverse();	
					}					
				}				
			}
		}
		return nextLinker;
	}
	
	void AnnotationLoops::dumpLoop(
		std::ostringstream& oss, 
		const Loop& aLoop) const
	{
		std::vector< Linker >::const_iterator it;
		for(it = aLoop.getLinkers().begin(); 
			it != aLoop.getLinkers().end(); 
			++ it)
		{
			oss << "{";
			if(0 < it->getResidues().size())
			{
				oss << it->getResidues().front();
				oss << "-";
				oss << it->getResidues().back();
			}
			else
			{
				oss << "no residues";
			}
			oss << "}, ";
		}
		oss << aLoop.describe();
	}	
	
	std::string AnnotationLoops::output() const
	{
		std::ostringstream oss;
		int i = 0;
		std::vector< Loop >::const_iterator it;
		for(it = mLoops.begin(); it != mLoops.end(); ++it)
		{
			oss << "Loop " << i << " : ";
			dumpLoop(oss, *it);
			oss << std::endl;
			++ i;
		}
		return oss.str();
	}
}
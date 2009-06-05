#include "AnnotationLoops.h"
#include "AnnotateModel.h"
#include "AnnotationLinkers.h"
#include <sstream>

namespace annotate
{
	AnnotationLoops::AnnotationLoops() 
	{
		addRequirement(AnnotationLinkers().provides());
	}
	
	AnnotationLoops::~AnnotationLoops() 
	{
		clear();
	}
	
	void AnnotationLoops::clear()
	{
		mLoops.clear();
	}
	
	const std::string AnnotationLoops::provides() const 
	{
		std::string strAnnotationName = "Loops";
		return strAnnotationName;
	}
	
	const std::vector< Loop >& AnnotationLoops::getLoops() const
	{
		return mLoops;
	}
	
	void AnnotationLoops::update(const AnnotateModel& aModel)
	{
		std::map<ResId, const Linker*> residueLinkerMap;
		residueLinkerMap = getResidueLinkerMap(aModel);
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
				if(potentialLoop.front().getStart().isValid() && potentialLoop.back().getEnd().isValid())
				{
					const Stem* startStem = &potentialLoop.front().getStart().getStem();
					const Stem* endStem = &potentialLoop.back().getEnd().getStem();
					if(startStem == endStem)
					{
						mLoops.push_back(Loop(potentialLoop));
					}
				}
			}
		}
	}
	
	std::map<mccore::ResId, const Linker*>
	AnnotationLoops::getResidueLinkerMap(const AnnotateModel& aModel) const
	{
		std::map<mccore::ResId, const Linker*> residueLinkerMap;
		const Annotation* pAnnotation = aModel.getAnnotation("Linkers");
		const AnnotationLinkers* pAnnotLinkers = NULL;
		pAnnotLinkers = dynamic_cast<const AnnotationLinkers*>(pAnnotation);
		if(NULL != pAnnotLinkers)
		{
			std::vector<Linker>::const_iterator it;
			for(it = pAnnotLinkers->getLinkers().begin();
				it != pAnnotLinkers->getLinkers().end(); 
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
	
	// TODO : This should be in the stem
	mccore::ResId AnnotationLoops::nextId(
		const Stem& aStem, 
		const StemConnection& aConnection) const
	{
		ResId id;
		switch(aConnection.getConnection())
		{
		case Stem::eFIRST_STRAND_FRONT_PAIR:
			id = aStem.basePairs().front().rResId;
			break;
		case Stem::eSECOND_STRAND_FRONT_PAIR:
			id = aStem.basePairs().front().fResId;
			break;
		case Stem::eFIRST_STRAND_BACK_PAIR:
			id = aStem.basePairs().back().rResId;
			break;
		case Stem::eSECOND_STRAND_BACK_PAIR:
			id = aStem.basePairs().back().fResId;
			break;
		default:
			// TODO : This should be an exception
			break;
		}
		
		return id;
	}
	
	Linker AnnotationLoops::nextLinker(
		const Linker& aLinker,
		const std::map<mccore::ResId, const Linker*>& aResidueLinkerMap) const
	{
		Linker nextLinker;
		const Stem* pNextStem = &aLinker.getEnd().getStem();
		if(NULL != pNextStem)
		{
			ResId resId = nextId(
				aLinker.getEnd().getStem(), 
				aLinker.getEnd());
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
			oss << it->getResidues().front();
			oss << "-";
			oss << it->getResidues().back();
			oss << "}, ";
		}
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
#include "AnnotationLoops.h"
#include "AnnotateModel.h"
#include "AnnotationLinkers.h"
#include <sstream>
#include <cassert>

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
	
	void AnnotationLoops::update(AnnotateModel& aModel)
	{
		const AnnotationLinkers* pAnnotLinkers = NULL;
		pAnnotLinkers = aModel.getAnnotation<AnnotationLinkers>();
		
		// Get the potential loops from the linkers
		std::list<Loop> potentialLoops;
		getIncompleteLoops(aModel, potentialLoops);
		
		while(!potentialLoops.empty())
		{
			bool bConnectionFound = false;
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
					bConnectionFound = true;
				}
				else
				{
					it ++;
				}
			}
			if(!bConnectionFound)
			{
				// Check loop integrity
				if(!workLoop.checkIntegrity())
				{
					mccore::gOut(0) << "Loop failed integrity check" << std::endl;
					std::vector<Linker>::const_iterator itLink;
					for(itLink = workLoop.linkers().begin(); 
						itLink != workLoop.linkers().end(); 
						++ itLink)
					{
						mccore::gOut(0) << "Linker : " << std::endl;
						mccore::gOut(0) << "residues : ";
						std::vector<mccore::ResId>::const_iterator itRes;
						for(itRes = itLink->residues().begin(); 
							itRes != itLink->residues().end();
							++ itRes)
						{
							if(itRes != itLink->residues().begin())
							{
								mccore::gOut(0) << "-";
							}
							mccore::gOut(0) << *itRes;
						}
						mccore::gOut(0) << std::endl;
					}
					// assert(false);
				}
				mLoops.push_back(workLoop);
			}
			else
			{
				// Loop connected, but could connect with more
				potentialLoops.push_back(workLoop);
			}
		}
	}
	
	void AnnotationLoops::getIncompleteLoops(
		const AnnotateModel& aModel, 
		std::list<Loop>& aLoops) const
	{
		const AnnotationLinkers* pALinkers = NULL;
		pALinkers = aModel.getAnnotation<AnnotationLinkers>();
		
		if(NULL != pALinkers)
		{
			std::vector<Linker>::const_iterator itLinker;
			for(itLinker = pALinkers->getLinkers().begin();
				itLinker != pALinkers->getLinkers().end(); 
				++ itLinker)
			{
				aLoops.push_back(Loop(*itLinker));
			}
		}		
	}
	
	bool AnnotationLoops::loopSelfComplete(const Loop& aLoop) const
	{
		bool bComplete = aLoop.closed() || aLoop.opened();
		return bComplete;
	}
	
	std::vector< Loop > AnnotationLoops::getLoops(
		const std::string& aDescription) const
	{
		std::vector< Loop > loops;
		std::vector< Loop >::const_iterator it;
		for(it = mLoops.begin(); it != mLoops.end(); ++ it)
		{
			if(0 == aDescription.compare(it->describe()))
			{
				loops.push_back(*it);
			}
		}
		return loops;
	}
	
	std::string AnnotationLoops::outputLoop(const Loop& aLoop) const
	{
		std::ostringstream oss;
		std::vector< Linker >::const_iterator it;
		for(it = aLoop.linkers().begin(); it != aLoop.linkers().end(); ++ it)
		{
			oss << "{";
			if(0 < it->residues().size())
			{
				oss << it->residues().front();
				oss << "-";
				oss << it->residues().back();
			}
			else
			{
				oss << "no residues";
			}
			oss << "}, ";
		}
		oss << aLoop.describe();
		return oss.str();
	}
	
	std::string AnnotationLoops::output() const
	{
		std::ostringstream oss;
		int i = 0;
		std::vector< Loop >::const_iterator it;
		for(it = mLoops.begin(); it != mLoops.end(); ++it)
		{
			oss << "Loop " << i << " : " << outputLoop(*it) << std::endl;
			++ i;
		}
		return oss.str();
	}
}
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
		loop_list potentialLoops = createLoopsFromLinkers(aModel);

		// Connect together the potential loops
		loop_list connectedLoops = findLoopsByConnectivity(potentialLoops);

		// Keep the final loops
		loop_list::const_iterator it;
		for(it = connectedLoops.begin(); it != connectedLoops.end(); ++ it)
		{
			mLoops.push_back(*it);
		}
	}

	std::list<Loop> AnnotationLoops::createLoopsFromLinkers(
				const AnnotateModel& aModel) const
	{
		std::list<Loop> loops;
		const AnnotationLinkers* pALinkers = NULL;
		pALinkers = aModel.getAnnotation<AnnotationLinkers>();

		if(NULL != pALinkers)
		{
			std::vector<Linker>::const_iterator itLinker;
			for(itLinker = pALinkers->getLinkers().begin();
				itLinker != pALinkers->getLinkers().end();
				++ itLinker)
			{
				loops.push_back(Loop(*itLinker));
			}
		}
		return loops;
	}

	std::list<Loop> AnnotationLoops::findLoopsByConnectivity(
			std::list<Loop>& aPotentials) const
	{
		std::list<Loop> loops;
		std::list<Loop> potentialLoops = aPotentials;

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
				// Loop is not connected to any other, this is a complete loop
				// considering only connectivity
				loops.push_back(workLoop);
			}
			else
			{
				// Loop connected, but could connect with more
				potentialLoops.push_back(workLoop);
			}
		}
		return loops;
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

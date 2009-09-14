#include "AnnotationLoops.h"
#include "AnnotateModel.h"
#include "AnnotationInteractions.h"
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
		addRequirement<AnnotationInteractions>();
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

		loop_list loops = findLoops(aModel, pAnnotLinkers->getLinkers());
		mLoops.insert(mLoops.end(), loops.begin(), loops.end());
	}

	std::list<Loop> AnnotationLoops::findLoops(
		AnnotateModel& aModel,
		const std::vector<AnnotationLinkers::linker_info>& aLinkers) const
	{
		const AnnotationInteractions* pAnnotInters = NULL;
		pAnnotInters = aModel.getAnnotation<AnnotationInteractions>();
		assert(NULL != pAnnotInters);

		std::list<Linker> looses;
		std::list<Linker> used;
		std::vector<AnnotationLinkers::linker_info>::const_iterator it;
		for(it = aLinkers.begin(); it != aLinkers.end(); ++ it)
		{
			if(it->second == NULL)
			{
				// Pseudo linker don't form open loops
				looses.push_back(it->first);
			}
		}

		std::list<Loop> loops;
		for(it = aLinkers.begin(); it != aLinkers.end(); ++ it)
		{
			if(it->second == NULL && NULL != it->first.start())
			{
				// Only start with real linkers and ignore open linkers for now
				std::list<Loop> loop = findLoop(*pAnnotInters, aLinkers, it);
				assert(loop.size() <= 1);
				if(loop.size() == 1)
				{
					loops.push_back(loop.front());
					used.insert(used.end(), loop.front().linkers().begin(), loop.front().linkers().end());
				}
			}
		}

		removeLooses(used, looses);
		std::list<Linker>::const_iterator itLoose;
		for(itLoose = looses.begin(); itLoose != looses.end(); ++ itLoose)
		{
			loops.push_back(Loop(*itLoose));
		}
		return loops;
	}

	std::list<Loop> AnnotationLoops::findLoop(
		const AnnotationInteractions& aInteractions,
		const std::vector<AnnotationLinkers::linker_info>& aLinkers,
		const std::vector<AnnotationLinkers::linker_info>::const_iterator& aIt) const
	{
		std::list<Loop> loops;
		std::vector<Linker> potentialLoop;
		std::vector<AnnotationLinkers::linker_info>::const_iterator it = aIt;
		potentialLoop.push_back(it->first);
		assert(NULL != it->first.start());
		const Stem* pStem = dynamic_cast<const Stem*>(it->first.start());
		if(NULL != pStem)
		{
			BasePair startPair = pStem->basePairs().back();
			if(startPair.fResId == it->first.residues().front())
			{
				++ it; // Start at next loop
				while(it != aLinkers.end() && potentialLoop.back().residues().back() != startPair.rResId)
				{
					if(!aInteractions.areContiguous(potentialLoop.back().residues().back(), it->first.residues().front())
							&& (potentialLoop.back().residues().back() != it->first.residues().front()))
					{
						break;
					}

					if(it->second != NULL)
					{
						if(it->second == pStem)
						{
							if(it->first.residues().front() == startPair.rResId)
							{
								// Completed a loop
								break;
							}
							else
							{
								// If this is false, we have a real knot !!!
								assert(false);
							}
						}
						else if(pStem->pseudoKnots(*it->second))
						{
							// Pseudo-knot, keep walking
							potentialLoop.push_back(it->first);
						}
						else
						{
							if(it->first.start() != NULL)
							{
								// Walk until we jumped to next strand
								it = advance(aLinkers, it, it->second->basePairs().front().rResId);
								if(it != aLinkers.end())
								{
									potentialLoop.push_back(it->first);
								}
								else
								{
									// Could not reach next strand, this is open
									break;
								}
							}
							else
							{
								// There was a cut between linkers
								break;
							}
						}
					}
					else
					{
						if(it->first.residues().front() == potentialLoop.back().residues().back())
						{
							potentialLoop.push_back(it->first);
						}
						else
						{
							// The new linker does not follow the previous one on the backbone
							break;
						}
					}
					++ it;
				}
				if(potentialLoop.back().residues().back() == startPair.rResId)
				{
					loops.push_back(Loop(potentialLoop));
					/*
					std::cout << "Loop completed : ";
					std::vector<Linker>::const_iterator itLinker;
					for(itLinker = potentialLoop.begin(); itLinker != potentialLoop.end(); ++ itLinker)
					{
						if(itLinker != potentialLoop.begin())
						{
							std::cout << ",";
						}
						std::cout << "{" << itLinker->residues().front() << "," << itLinker->residues().back() << "}";
					}
					std::cout << std::endl;*/
				}
			}
		}
		return loops;
	}

	std::vector<AnnotationLinkers::linker_info>::const_iterator
	AnnotationLoops::advance(
		const std::vector<AnnotationLinkers::linker_info>& aLinkers,
		const std::vector<AnnotationLinkers::linker_info>::const_iterator& aIt,
		const LabeledResId& aResId) const
	{
		std::vector<AnnotationLinkers::linker_info>::const_iterator it = aIt;
		while(it != aLinkers.end() && it->first.residues().front() != aResId)
		{
			++ it;
		}
		return it;
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

	void AnnotationLoops::removeLooses(
		const std::list<Linker>& aUsed,
		std::list<Linker>& aLooses) const
	{
		std::list<Linker>::iterator itLoose = aLooses.begin();
		std::list<Linker>::const_iterator itUsed = aUsed.begin();
		// TODO : This should be O(N)
		for(itUsed = aUsed.begin(); itUsed != aUsed.end(); ++ itUsed)
		{
			for(itLoose = aLooses.begin(); itLoose != aLooses.end(); ++ itLoose)
			{
				if(*itLoose == *itUsed)
				{
					aLooses.erase(itLoose);
					break;
				}
			}
		}
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

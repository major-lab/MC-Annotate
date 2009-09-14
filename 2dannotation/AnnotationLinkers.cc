#include "AnnotationLinkers.h"
#include "AnnotateModel.h"
#include "AnnotationInteractions.h"
#include "AnnotationStems.h"
#include <cassert>
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationLinkers::mstrAnnotationName = "Linkers";

	// Methods
	AnnotationLinkers::AnnotationLinkers()
	{
		addRequirement<AnnotationInteractions>();
		addRequirement<AnnotationStems>();
	}

	AnnotationLinkers::~AnnotationLinkers()
	{
		clear();
	}

	void AnnotationLinkers::clear()
	{
		mResidueInfos.clear();
		mLinkers.clear();
	}

	void AnnotationLinkers::update(AnnotateModel& aModel)
	{
		clear();

		computeResidueInfos(aModel);

		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>();

		if(NULL != pAnnotStems && 0 < pAnnotStems->getStems().size())
		{
			updateLinkers();
		}
	}

	void AnnotationLinkers::updateLinkers()
	{
		std::vector< std::vector<stResidueInfo> >::const_iterator it;
		for(it = mResidueInfos.begin(); it != mResidueInfos.end(); ++ it)
		{
			updateChainLinkers(*it);
		}

		// Chain the pseudo linkers to other linkers
		connectLinkersToLinkers();
	}

	std::list<std::list<unsigned int> > AnnotationLinkers::getLinkerRanges(
		const std::vector<stResidueInfo>& chainInfo) const
	{
		const Stem* pPrevStem = NULL;
		std::list<std::list<unsigned int> > ranges;
		std::list<unsigned int> residues;
		for(unsigned int uiIndex = 0; uiIndex < chainInfo.size(); ++ uiIndex)
		{
			if(pPrevStem == chainInfo[uiIndex].pStem)
			{
				residues.push_back(uiIndex);
			}
			else if(0 < residues.size())
			{
				assert(0 < residues.size());
				pPrevStem = chainInfo[uiIndex].pStem;
				ranges.push_back(residues);
				residues.clear();
				residues.push_back(uiIndex);
			}
			else
			{
				// Starting directly on a stem
				assert(0 == uiIndex);
				residues.push_back(uiIndex);
				pPrevStem = chainInfo[uiIndex].pStem;
			}
		}
		if(0 < residues.size())
		{
			ranges.push_back(residues);
		}
		return ranges;
	}

	std::list<AnnotationLinkers::linker_info> AnnotationLinkers::getLinker(
		const std::vector<stResidueInfo>& aChainInfo,
		const std::list<unsigned int>& aRange)
	{
		assert(0 < aRange.size());
		std::list<AnnotationLinkers::linker_info> linkers;
		const Stem* pStem = aChainInfo[aRange.front()].pStem;
		const Stem* pPrev = NULL;
		const Stem* pNext = NULL;
		std::pair<Linker, const Stem*> linkerInfo;
		std::vector<LabeledResId> linkerResidues;
		std::list<unsigned int>::const_iterator itRes;
		for(itRes = aRange.begin(); itRes != aRange.end(); ++ itRes)
		{
			linkerResidues.push_back(aChainInfo[*itRes].resId);
		}
		if(NULL == pStem)
		{
			if(0 < aRange.front())
			{
				linkerResidues.insert(linkerResidues.begin(), aChainInfo[aRange.front() - 1].resId);
				pPrev = aChainInfo[aRange.front() - 1].pStem;
			}
			if((aRange.back() + 1) < aChainInfo.size())
			{
				linkerResidues.push_back(aChainInfo[aRange.back() + 1].resId);
				pNext = aChainInfo[aRange.back() + 1].pStem;
			}
		}
		if(2 <= linkerResidues.size() && ((NULL != pNext || NULL != pPrev) || NULL != pStem))
		{
			linkerInfo.first = createLinker(linkerResidues, pPrev, pNext);
			linkerInfo.second = pStem;
			linkers.push_back(linkerInfo);
		}
		else
		{
			// we've decided to drop this one
			if(0 <= linkerResidues.size())
			{
				mUnconnected.insert(linkerResidues.front());
				mUnconnected.insert(linkerResidues.back());
			}
		}
		return linkers;
	}

	bool AnnotationLinkers::shouldConnect(
		const linker_info& aFirst,
		const linker_info& aSecond) const
	{
		bool bShouldConnect = false;

		LabeledResId firstResId = aFirst.first.residues().back();
		LabeledResId secondResId = aSecond.first.residues().front();
		bShouldConnect = (mUnconnected.end() == mUnconnected.find(firstResId)
			&& mUnconnected.end() == mUnconnected.find(secondResId));
		return bShouldConnect;
	}

	void AnnotationLinkers::connectLinkersToLinkers()
	{
		std::vector<linker_info>::iterator it;
		std::vector<linker_info>::iterator itPrev = mLinkers.end();
		std::vector<linker_info>::iterator itNext;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++it)
		{
			itNext = it;
			++ itNext;
			if(NULL == it->first.start() && itPrev != mLinkers.end())
			{
				if(shouldConnect(*itPrev, *it))
				{
					it->first.start(&itPrev->first);
				}
			}
			if(NULL == it->first.end() && itNext != mLinkers.end())
			{
				if(shouldConnect(*it, *itNext))
				{
					it->first.end(&itNext->first);
				}
			}
			itPrev = it;
		}
	}

	void AnnotationLinkers::connectPseudoLinkers()
	{
		// Add the missing links
		std::vector<linker_info>::iterator it = mLinkers.begin();
		std::vector<linker_info>::iterator itNext;
		while(it != mLinkers.end())
		{
			itNext = it;
			itNext ++;
			if(mLinkers.end() != itNext)
			{
				if(it->second != NULL && itNext->second != NULL)
				{
					// Two unconnected pseudo-linkers
					if(shouldConnect(*it, *itNext))
					{
						linker_info newLinker;
						std::vector<LabeledResId> linkerResidues;
						linkerResidues.push_back(it->first.residues().back());
						linkerResidues.push_back(itNext->first.residues().front());
						newLinker.first = createLinker(linkerResidues, it->second, itNext->second);
						newLinker.second = NULL;
						it = mLinkers.insert(itNext, newLinker);
					}
				}
			}
			it ++;
		}
	}

	void AnnotationLinkers::updateChainLinkers(
		const std::vector<stResidueInfo>& chainInfo)
	{
		std::vector<linker_info> linkers;
		std::list<std::list<unsigned int> > ranges;

		ranges = getLinkerRanges(chainInfo);

		std::list<std::list<unsigned int> >::const_iterator it;
		for(it = ranges.begin(); it != ranges.end(); ++ it)
		{
			std::list<linker_info> linker = getLinker(chainInfo, *it);
			if(1 == linker.size())
			{
				mLinkers.push_back(linker.front());
			}
		}

		// Make the connections with the pseudo-linkers
		connectPseudoLinkers();
	}

  	Linker AnnotationLinkers::createLinker(
  		const std::vector<LabeledResId>& aResidues,
  		const SecondaryStructure* apStartStruct,
  		const SecondaryStructure* apEndStruct) const
  	{
  		assert(2 <= aResidues.size() || (NULL == apStartStruct && NULL == apEndStruct));
		Linker linker(aResidues, apStartStruct, apEndStruct);
		if(2 <= aResidues.size())
		{
			linker.order();
		}
		return linker;
  	}

  	void AnnotationLinkers::computeResidueInfos(const AnnotateModel& aModel)
	{
  		const AnnotationInteractions* pAnnotInters = aModel.getAnnotation<AnnotationInteractions>();
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>();
		assert(NULL != pAnnotInters);
		assert(NULL != pAnnotStems);

		std::vector<stResidueInfo> chainInfo;
		mccore::GraphModel::const_iterator it = aModel.begin();
		char cPrevChain = it->getResId().getChainId();
		for(;it != aModel.end(); ++ it)
		{
			if(cPrevChain != it->getResId().getChainId())
			{
				if(0 < chainInfo.size())
				{
					mUnconnected.insert(chainInfo.back().resId);
					mUnconnected.insert(it->getResId());
					mResidueInfos.push_back(chainInfo);
					chainInfo.clear();
				}
				cPrevChain = it->getResId().getChainId();
			} else if(0 < chainInfo.size() && !pAnnotInters->areContiguous(chainInfo.back().resId, it->getResId()))
			{
				mUnconnected.insert(chainInfo.back().resId);
				mUnconnected.insert(it->getResId());
				mResidueInfos.push_back(chainInfo);
				chainInfo.clear();
			}

			stResidueInfo resInfo;
			resInfo.resId = (*it).getResId();
			resInfo.pStem = NULL;
			if(NULL != pAnnotStems)
			{
				std::vector<Stem>::const_iterator stemIt;
				for(stemIt = pAnnotStems->getStems().begin();
					stemIt != pAnnotStems->getStems().end();
					++stemIt)
				{
					if((*stemIt).contains((*it).getResId()))
					{
						resInfo.pStem = &(*stemIt);
					}
				}
			}
			chainInfo.push_back(resInfo);
		}
		if(0 < chainInfo.size())
		{
			mResidueInfos.push_back(chainInfo);
		}
	}

	std::string AnnotationLinkers::outputLinker(const Linker& aLinker) const
	{
		std::ostringstream oss;
		if(0 < aLinker.residues().size())
		{
			oss << aLinker.residues().front();
			oss << "-";
			oss << aLinker.residues().back();
		}else
		{
			oss << "no residues";
		}
		return oss.str();
	}

	std::string AnnotationLinkers::output() const
	{
		std::ostringstream oss;
		int i = 0;
		std::vector<linker_info>::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++it)
		{
			if(NULL == it->second)
			{
				oss << "Linker " << i << " : " << outputLinker(it->first) << std::endl;
				++ i;
			}
		}
		i = 0;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++ it)
		{
			if(NULL != it->second)
			{
				oss << "Pseudo linker " << i << " : " << outputLinker(it->first) << std::endl;
				++ i;
			}
		}
		return oss.str();
	}
}

#include "AnnotationLinkers.h"
#include "AnnotateModel.h"
#include "AnnotationStems.h"
#include <sstream>

namespace annotate
{
	// Static members
	std::string AnnotationLinkers::mstrAnnotationName = "Linkers";
	
	// Methods
	AnnotationLinkers::AnnotationLinkers()
	{
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
		
	void AnnotationLinkers::update(const AnnotateModel& aModel)
	{
		clear();
		
		computeResidueInfos(aModel);
		
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>();
		
		if(NULL != pAnnotStems)
		{
			std::set<Linker> linkerSet;
			if(0 < pAnnotStems->getStems().size())
			{
				updateLinkers(linkerSet);
			}
			else if(0 < aModel.size())
			{
				std::vector< std::vector<stResidueInfo> >::const_iterator it;
				for(it = mResidueInfos.begin(); it != mResidueInfos.end(); ++ it)
				{
					allResiduesLinker(*it, linkerSet);
				}
			}
			
			std::set<Linker>::const_iterator itLinker = linkerSet.begin();
			for(;itLinker != linkerSet.end(); ++itLinker)
			{
				mLinkers.push_back(*itLinker);
			}
		}
	}
	
	void AnnotationLinkers::updateLinkers(std::set<Linker>& linkers) const
	{
		std::vector< std::vector<stResidueInfo> >::const_iterator it;
		for(it = mResidueInfos.begin(); it != mResidueInfos.end(); ++ it)
		{
			updateChainLinkers(*it, linkers);
		}
	}
	
	void AnnotationLinkers::updateChainLinkers(
		const std::vector<stResidueInfo>& chainInfo, 
		std::set<Linker>& linkers) const
	{
		const Stem* pPrevStem = NULL;
		std::vector<stResidueInfo>::const_iterator it;
		std::vector<stResidueInfo>::const_iterator itPrevStem = chainInfo.end();
		std::vector<mccore::ResId> linkerResidues;
		for(it = chainInfo.begin(); it != chainInfo.end(); ++ it)
		{
			if(NULL == pPrevStem && NULL == it->pStem)
			{
				// First linker
				linkerResidues.push_back(it->resId);
			}
			else if(NULL == pPrevStem && NULL != it->pStem)
			{
				// First linker completed
				if(0 < linkerResidues.size())
				{
					Linker linker = createLinker(
						linkerResidues, 
						NULL, mccore::ResId(), 
						it->pStem, it->resId);
					linkers.insert(linker);
					linkerResidues.clear();  // Starting a new linker
				}
				pPrevStem = it->pStem;
				itPrevStem = it;
			}
			else if(NULL != pPrevStem && it->pStem == pPrevStem)
			{
				// Walking a stem
				if(0 < linkerResidues.size())
				{
					// This is the end of an hairpin loop
					Linker linker = createLinker(
						linkerResidues, 
						it->pStem, itPrevStem->resId, 
						it->pStem, it->resId);
					linkers.insert(linker);
					linkerResidues.clear();  // Starting a new linker
				}
				itPrevStem = it;	
			}
			else if(NULL != pPrevStem && it->pStem == NULL)
			{
				// Between stems
				linkerResidues.push_back(it->resId);
			}
			else if(NULL != pPrevStem && it->pStem != NULL)
			{
				if(pPrevStem->pseudoKnots(*it->pStem))
				{
					// Keep walking
					linkerResidues.push_back(it->resId);
				}
				else
				{
					// Connection between stems done
					Linker linker = createLinker(
						linkerResidues, 
						pPrevStem, itPrevStem->resId, 
						it->pStem, it->resId);
					linkers.insert(linker);
					linkerResidues.clear();  // Starting a new linker
					pPrevStem = it->pStem;
				}
			}
		}
		if(0 < linkerResidues.size())
		{
			mccore::ResId resId;
			const Stem* pStem = NULL;
			if(itPrevStem != chainInfo.end())
			{
				resId = itPrevStem->resId;
				pStem = itPrevStem->pStem;
			}
			Linker linker = createLinker(
				linkerResidues, 
				pStem, resId, 
				NULL, mccore::ResId());
			linkers.insert(linker);
			linkerResidues.clear();  // Starting a new linker
		}
	}
	
  	void AnnotationLinkers::allResiduesLinker(
  		const std::vector<stResidueInfo>& chainInfo, 
		std::set<Linker>& linkers) const
  	{
  		std::vector<mccore::ResId> linkerResidues;
  		std::vector<stResidueInfo>::const_iterator it;
  		for(it = chainInfo.begin(); it != chainInfo.end(); ++ it)
  		{
  			linkerResidues.push_back(it->resId);
  		}
  		Linker linker = Linker(linkerResidues, StemConnection(), StemConnection());
	  	linker.order();
	  	linkers.insert(linker);	  	
  	}
  	
  	Linker AnnotationLinkers::createLinker(
  		const std::vector<mccore::ResId>& aResidues, 
  		const Stem* apStem1, 
  		const mccore::ResId& aResId1, 
  		const Stem* apStem2, 
  		const mccore::ResId& aResId2) const
  	{
  		StemConnection startConnect;
  		StemConnection endConnect;
  		if(NULL != apStem1)
  		{
  			Stem::enConnection eStartConnect = apStem1->getConnection(aResId1);
  			startConnect = StemConnection(*apStem1, eStartConnect); 
  		}
  		if(NULL != apStem2)
  		{
  			Stem::enConnection eEndConnect = apStem2->getConnection(aResId2);
  			endConnect = StemConnection(*apStem2, eEndConnect); 
  		}
  		
		Linker linker(aResidues, startConnect, endConnect);
		linker.order();
		return linker;
  	}
  	
  	void AnnotationLinkers::computeResidueInfos(const AnnotateModel& aModel)
	{
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>();
		
		std::vector<stResidueInfo> chainInfo;
		GraphModel::const_iterator it = aModel.begin();
		char cPrevChain = it->getResId().getChainId();
		for(;it != aModel.end(); ++ it)
		{
			if(cPrevChain != it->getResId().getChainId())
			{
				mResidueInfos.push_back(chainInfo);
				chainInfo.clear();
				cPrevChain = it->getResId().getChainId();
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
		if(0 < aLinker.getResidues().size())
		{
			oss << aLinker.getResidues().front();
			oss << "-"; 
			oss << aLinker.getResidues().back();
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
		std::vector<Linker>::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++it)
		{
			oss << "Linker " << i << " : " << outputLinker(*it) << std::endl;
			/*
			if(0 < it->getResidues().size())
			{
				oss << it->getResidues().front();
				oss << "-"; 
				oss << it->getResidues().back();
			}else
			{
				oss << "no residues";
			}
			oss << std::endl;*/
			++ i;
		}
		return oss.str();
	}
		
	const std::vector< Linker >& AnnotationLinkers::getLinkers() const
	{
		return mLinkers;
	}
}
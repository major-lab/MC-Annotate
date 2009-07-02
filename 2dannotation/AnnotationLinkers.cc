#include "AnnotationLinkers.h"
#include "AnnotateModel.h"
#include "AnnotationStems.h"
#include <sstream>

namespace annotate
{
	AnnotationLinkers::AnnotationLinkers()
	{
		addRequirement(AnnotationStems().provides());
	}
	
	AnnotationLinkers::~AnnotationLinkers()
	{
		clear();
	}
	
	const std::string AnnotationLinkers::provides() const
	{
		std::string strAnnotationName = "Linkers";
		return strAnnotationName;
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
		
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>(AnnotationStems().provides());
		
		if(NULL != pAnnotStems)
		{
			std::set<Linker> linkerSet;
			if(0 < pAnnotStems->getStems().size())
			{		
			  	std::vector< Stem >::const_iterator it;
				for(it = pAnnotStems->getStems().begin(); 
					it != pAnnotStems->getStems().end(); 
					++it)
				{
					findLinker((&*it), Stem::eFIRST_STRAND_FRONT_PAIR, linkerSet);
					findLinker((&*it), Stem::eFIRST_STRAND_BACK_PAIR, linkerSet);
					findLinker((&*it), Stem::eSECOND_STRAND_FRONT_PAIR, linkerSet);
					findLinker((&*it), Stem::eSECOND_STRAND_BACK_PAIR, linkerSet);		
				}
			}
			else if(0 < aModel.size())
			{
				linkerSet.insert(allResiduesLinker());
			}
			
			std::set<Linker>::const_iterator itLinker = linkerSet.begin();
			for(;itLinker != linkerSet.end(); ++itLinker)
			{
				mLinkers.push_back(*itLinker);
			}
		}
	}
	
	Linker AnnotationLinkers::findLinker(const StemConnection& aConnection) const
	{
	  	Linker linker;
	  	
	  	// Find the possible strand from this end
	  	int iDirection = aConnection.getDirection();
	  	ResId resId = aConnection.getResidue();
	
	  	// Find the residue information
	  	std::vector<stResidueInfo>::const_iterator it = findResidueInfo(resId);
	  	
	  	std::vector<mccore::ResId> linkerResidues;
	  	
	  	// Set the start of the linker
	  	StemConnection linkerStart = aConnection;
	 	StemConnection linkerEnd;
	  	
	  	// Walk the chain until we find a stem not pseudo-knotted
	  	if(0 < iDirection)
	  	{
		  	std::advance(it, iDirection);
		  	while(it != mResidueInfos.end())
		  	{
		  		const Stem* pStem = (*it).pStem;
		  		if(NULL != pStem && !linkerStart.getStem().pseudoKnots(*pStem))
		  		{
		  			// Stem found
		  			Stem::enConnection eConnect = pStem->getConnection((*it).resId);
		  			linkerEnd = StemConnection(*pStem, eConnect);
		  			it = mResidueInfos.end();
		  		}
		  		else
		  		{
		  			linkerResidues.push_back((*it).resId);
		  			std::advance(it, iDirection);
		  		}
		  	}
	  	}
	  	linker = Linker(linkerResidues, linkerStart, linkerEnd);
	  	linker.order();
	  	
	  	return linker;
	}
  
  	void AnnotationLinkers::findLinker(
  		const Stem* apStem, 
  		const Stem::enConnection& aeConnect,
  		std::set<Linker>& outLinkerSet) const
	{
	  	StemConnection connect(*apStem, aeConnect);
	  		
		Linker linker = findLinker(connect);
		
		if(!linker.isEmpty()) 
		{
			outLinkerSet.insert(linker); 
		}
  	}
  	
  	Linker AnnotationLinkers::allResiduesLinker() const
  	{
  		std::vector<mccore::ResId> linkerResidues;
  		std::vector<stResidueInfo>::const_iterator it;
  		for(it = mResidueInfos.begin(); it != mResidueInfos.end(); ++ it)
  		{
  			linkerResidues.push_back(it->resId);
  		}
  		Linker linker = Linker(linkerResidues, StemConnection(), StemConnection());
	  	linker.order();
	  	return linker;
  	}
  	
  	void AnnotationLinkers::computeResidueInfos(const AnnotateModel& aModel)
	{
		int i = 0;
		const AnnotationStems* pAnnotStems = aModel.getAnnotation<AnnotationStems>(AnnotationStems().provides());
		
		if(NULL != pAnnotStems)
		{
			mResidueInfos.resize(aModel.size());
			GraphModel::const_iterator it = aModel.begin();
			for(;it != aModel.end(); ++ it)
			{
				mResidueInfos[i].resId = (*it).getResId();
				mResidueInfos[i].pStem = NULL;
				std::vector<Stem>::const_iterator stemIt;
				for(stemIt = pAnnotStems->getStems().begin(); 
					stemIt != pAnnotStems->getStems().end(); 
					++stemIt)
				{
					if((*stemIt).contains((*it).getResId()))
					{
						if(NULL != mResidueInfos[i].pStem)
						{
							gOut(1) << "Residue " << (*it).getResId();
							gOut(1) << " associated with more than one stem";
							gOut(1) << std::endl;
						}
						mResidueInfos[i].pStem = &(*stemIt);
					}
				}
				++ i;
			}
		}	
	}
	
	std::vector<AnnotationLinkers::stResidueInfo>::const_iterator 
	AnnotationLinkers::findResidueInfo(const mccore::ResId& aResId) const
	{
		std::vector<AnnotationLinkers::stResidueInfo>::const_iterator it;
		for(it = mResidueInfos.begin(); it != mResidueInfos.end(); ++it)
		{
			if(aResId == (*it).resId)
			{
				break;
			}
		}
		return it;
	}
	
	std::string AnnotationLinkers::output() const
	{
		std::ostringstream oss;
		int i = 0;
		std::vector<Linker>::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++it)
		{
			oss << "Linker " << i << " : ";
			if(0 < it->getResidues().size())
			{
				oss << it->getResidues().front();
				oss << "-"; 
				oss << it->getResidues().back();
			}else
			{
				oss << "no residues";
			}
			oss << std::endl;
			++ i;
		}
		return oss.str();
	}
		
	const std::vector< Linker >& AnnotationLinkers::getLinkers() const
	{
		return mLinkers;
	}
}
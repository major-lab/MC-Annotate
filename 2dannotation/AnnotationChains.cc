#include "AnnotationChains.h"
#include "AnnotateModel.h"
#include <sstream>

namespace annotate
{
	std::string AnnotationChains::mstrAnnotationName = "Chains";
	
	AnnotationChains::AnnotationChains()
	{
	}
	
	AnnotationChains::~AnnotationChains()
	{
		clear();
	}
	
	void AnnotationChains::clear()
	{
		chain_map::iterator it;
		for(it = mChains.begin(); it != mChains.end(); ++ it)
		{
			it->second.clear();
		}
		mChains.clear();
	}
		
	void AnnotationChains::update(const AnnotateModel& aModel)
	{
		clear();
		
		chain_content chain;
		GraphModel::const_iterator it = aModel.begin();
		char cPrevChain = it->getResId().getChainId();
		for(;it != aModel.end(); ++ it)
		{
			if(cPrevChain != it->getResId().getChainId())
			{
				mChains.insert(chain_entry(cPrevChain, chain));
				chain.clear();
				cPrevChain = it->getResId().getChainId();
			}
			chain.push_back(it->getResId());
		}
		if(0 < chain.size())
		{
			mChains.insert(chain_entry(cPrevChain, chain));
		}
	}
	
	std::string AnnotationChains::output() const
	{
		ostringstream oss;
		
		chain_map::const_iterator itChain;
		for(itChain = mChains.begin(); mChains.end() != itChain; ++itChain)
	    {
	    	oss << "Chain " << itChain->first << " : ";
	    	chain_content::const_iterator itRes;
	    	for(itRes = itChain->second.begin(); 
	    		itRes != itChain->second.end(); 
	    		++ itRes)
	    	{
	    		if(itRes != itChain->second.begin())
	    		{
	    			oss << "-";
	    		}
	    		oss << *itRes;
	    	}
	    	oss << std::endl;
	    }
	  	
		return oss.str();		
	}
}
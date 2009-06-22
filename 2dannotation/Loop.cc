#include "Loop.h"

namespace annotate
{
	Loop::Loop() 
	{ 
	}
	
	Loop::Loop(const std::vector< Linker >& aLinkers)
	{
		mLinkers = aLinkers;
	}
	
	Loop::~Loop() 
	{
		clear();
	}
	
	const std::vector< Linker >& Loop::getLinkers() const
	{
		return mLinkers;
	}
	
	void Loop::clear()
	{
		mLinkers.clear();
	}
	
	void Loop::reverse()
	{
		std::vector< Linker >::iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++ it)
		{
			it->reverse();
		}
	}
	
	void Loop::append(const Loop& aLoop)
	{
		// Check that this can be added
		// TODO : Add exception
		
		mLinkers.insert(
			mLinkers.end(), 
			aLoop.mLinkers.begin(), 
			aLoop.mLinkers.end());
	}
	
	bool Loop::operator ==(const Loop& other) const
	{
		bool bEqual = false;
		
		if(mLinkers.size() == other.mLinkers.size())
		{
			bEqual = std::equal(
				mLinkers.begin(), 
				mLinkers.end(), 
				other.mLinkers.begin());
		}
		return bEqual;
	}
	
	bool Loop::contains(const mccore::ResId& aResId) const
	{
		bool bContains = false;
		std::vector< Linker >::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end() && !bContains; ++ it)
		{
			bContains = it->contains(aResId);
		}
		return bContains;
	}
	
	bool Loop::isAdjacent(const SecondaryStructure& aStruct) const
	{
		bool bAdjacent = false;
		
		const Loop* pLoop = dynamic_cast<const Loop*>(&aStruct);
		if(NULL != pLoop && operator == (*pLoop))
		{
			bAdjacent = true;
		}
		else
		{
						
			const Stem* pStem = dynamic_cast<const Stem*>(&aStruct);
			if(NULL != pStem)
			{
				// This is a stem
				std::vector<Linker>::const_iterator it = mLinkers.begin();
				for(; it != mLinkers.end() && !bAdjacent; ++ it)
				{
					bAdjacent = it->isAdjacent(*pStem);					
				}
			}
		}
		return bAdjacent;
	}
}
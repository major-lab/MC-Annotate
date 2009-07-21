#include "Loop.h"

namespace annotate
{
	Loop::Loop() 
	{ 
	}
	
	Loop::Loop(const Linker& aLinker)
	{
		mLinkers.push_back(aLinker);
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
		if(mLinkers.back().connects(aLoop.mLinkers.front()))
		{
			// Basic case, appends the linkers at the end of this one
			append_back(aLoop);
		}
		else if(mLinkers.back().connects(aLoop.mLinkers.back()))
		{
			// Flip the other linker before inserting it at the end of this one
			append_back_reverse(aLoop);			
		}
		else if(mLinkers.front().connects(aLoop.mLinkers.front()))
		{
			// Flip the other linker before inserting it at the beginning of this one
			append_front_reverse(aLoop);			
		}
		else if(mLinkers.front().connects(aLoop.mLinkers.back()))
		{
			// Insert the other loop in front of this one
			append_front(aLoop);
		}
		else
		{
			// TODO : Throw exception, this can't be added
			exit(0);
		}
		
		// Insure that the loop is still in orders
		order();
	}
	
	void Loop::append_back(const Loop& aLoop)
	{
		// Basic case, appends the linkers at the end of this one
		mLinkers.insert(
			mLinkers.end(), 
			aLoop.mLinkers.begin(), 
			aLoop.mLinkers.end());
	}
	
	void Loop::append_back_reverse(const Loop& aLoop)
	{
		// Flip the other linker before inserting it at the end of this one
			Loop workLoop = aLoop;
			workLoop.reverse();
			mLinkers.insert(
				mLinkers.end(), 
				workLoop.mLinkers.begin(), 
				workLoop.mLinkers.end());
	}
	
	void Loop::append_front(const Loop& aLoop)
	{
		// Insert the other loop in front of this one
		mLinkers.insert(
			mLinkers.begin(), 
			aLoop.mLinkers.begin(), 
			aLoop.mLinkers.end());
	}
	
	void Loop::append_front_reverse(const Loop& aLoop)
	{
		// Insert the other loop in front of this one
		Loop workLoop = aLoop;
		workLoop.reverse();
		mLinkers.insert(
			mLinkers.begin(), 
			workLoop.mLinkers.begin(), 
			workLoop.mLinkers.end());;
	}
	
	void Loop::order()
	{
		if(0 < mLinkers.size())
		{
			if(1 == mLinkers.size())
			{
				mLinkers.front().order();
			}
			else if(mLinkers.back() <  mLinkers.front())
			{
				reverse();
			}
		}
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
	
	std::string Loop::describe() const
	{
		std::string strDescription;
		// qualify the loop
		const StemConnection* pStart = &getLinkers().front().getStart();
		const StemConnection* pEnd = &getLinkers().back().getEnd();
		if(	pStart->isValid() 
			&& pEnd->isValid() 
			&& (&pStart->getStem() == &pEnd->getStem()))
		{
			int iSize = getLinkers().size();
			switch(iSize)
			{
			case 1:
				strDescription = "hairpin";
				break;
			case 2: 
				strDescription = "internal";
				break;
			default:
				strDescription = "multibranch";
			}
		}
		else
		{
			strDescription = "open";
		}
		return strDescription;
	}
	
	bool Loop::connects(const Loop& aLoop) const
	{
		bool bConnects = false;
		
		if(mLinkers.front().connects(aLoop.mLinkers.front()))
		{
			bConnects = true;
		}
		else if(mLinkers.front().connects(aLoop.mLinkers.back()))
		{
			bConnects = true;
		}
		else if(mLinkers.back().connects(aLoop.mLinkers.front()))
		{
			bConnects = true;
		}
		else if(mLinkers.back().connects(aLoop.mLinkers.back()))
		{
			bConnects = true;
		}
		return bConnects;		
	}
	
	bool Loop::complete() const
	{
		bool bComplete = false;
		
		if(0 < mLinkers.size())
		{
			if(mLinkers.front().connects(mLinkers.back()))
			{
				// Complete loop
				bComplete = true;
			}
			else
			{
				// Check for open loop completeness
				bool bStartConnected = mLinkers.front().getStart().isValid();
				bool bEndConnected = mLinkers.back().getEnd().isValid();
				bComplete = (!bStartConnected && !bEndConnected);
			}
		}
		return bComplete;		
	}
}
#include "Loop.h"

#include "BaseLink.h"
#include <cassert>

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

	bool Loop::isSame(const SecondaryStructure& aStruct) const
	{
		bool bSame = false;
		const Loop* pLoop = dynamic_cast<const Loop*>(&aStruct);
		if(NULL != pLoop && operator == (*pLoop))
		{
			// This is a loop and it has the same value as this one
			bSame = true;
		}
		return bSame;
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
		if(closed())
		{
			int iSize = mLinkers.size();
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

	bool Loop::closed() const
	{
		bool bClosed = false;

		if(0 < mLinkers.size())
		{
			bClosed = mLinkers.front().connects(mLinkers.back());
		}

		return bClosed;
	}

	bool Loop::complete() const
	{
		bool bComplete = false;
		if(closed())
		{
			// Ends are connected together
			bComplete = true;
		}
		else if(0 < mLinkers.size())
		{
			// ... or not connected at all
			bool bStartConnected = mLinkers.front().getStart().isValid();
			bool bEndConnected = mLinkers.back().getEnd().isValid();
			bComplete = (!bStartConnected && !bEndConnected);
		}
		return bComplete;
	}

	std::set<BaseInteraction> Loop::getBaseInteractions() const
	{
		std::set<BaseInteraction> interactions;
		std::vector<Linker>::const_iterator itLinker;
		for(itLinker = mLinkers.begin(); itLinker != mLinkers.end(); ++itLinker)
		{
			std::set<BaseInteraction> linkerInteractions;
			linkerInteractions = itLinker->getBaseInteractions();
			interactions.insert(linkerInteractions.begin(), linkerInteractions.end());
			assert(itLinker->getStart().isValid() || itLinker->getEnd().isValid());
			if(itLinker->getStart().isValid())
			{
				BasePair basePair = itLinker->getStart().getPair();
				BaseInteraction interact(
					basePair.first, basePair.fResId,
					basePair.second, basePair.rResId);
				interactions.insert(interact);
			}
			if(itLinker->getEnd().isValid())
			{
				BasePair basePair = itLinker->getEnd().getPair();
				BaseInteraction interact(
					basePair.first, basePair.fResId,
					basePair.second, basePair.rResId);
				interactions.insert(interact);
			}
		}
		return interactions;
	}

	void Loop::getLinkerInteractions(
		const Linker& aLinker,
		std::set<BaseLink>& aInteractions) const
	{
		std::set<BaseInteraction> links;
		links = aLinker.getBaseInteractions();
		std::set<BaseInteraction>::const_iterator it;
		for(it = links.begin(); it != links.end(); ++ it)
		{
			BaseLink link(it->first, it->fResId, it->second, it->rResId);
			aInteractions.insert(link);
		}
	}

	std::pair<std::set<BaseLink>, std::set<BasePair> > Loop::getInteractions() const
	{
		std::pair<std::set<BaseLink>, std::set<BasePair> > interactions;
		std::vector<Linker>::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++it)
		{
			getLinkerInteractions(*it, interactions.first);
			assert(it->getStart().isValid() || it->getEnd().isValid());
			if(it->getStart().isValid())
			{
				BasePair basePair = it->getStart().getPair();
				interactions.second.insert(basePair);
			}
			if(it->getEnd().isValid())
			{
				BasePair basePair = it->getEnd().getPair();
				interactions.second.insert(basePair);
			}
		}
		return interactions;
	}

	std::set<mccore::ResId> Loop::getResIds() const
	{
		std::set<mccore::ResId> resids;
		std::vector<Linker>::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++ it)
		{
			resids.insert(it->residues().begin(), it->residues().end());
			if(it->getStart().isValid())
			{
				BasePair basePair = it->getStart().getPair();
				resids.insert(basePair.fResId);
				resids.insert(basePair.rResId);
			}
			if(it->getEnd().isValid())
			{
				BasePair basePair = it->getEnd().getPair();
				resids.insert(basePair.fResId);
				resids.insert(basePair.rResId);
			}
		}
		return resids;
	}

	bool Loop::checkIntegrity() const
	{
		bool bIntegrity = true;
		std::vector<Linker>::const_iterator itLinker;
		for(itLinker = mLinkers.begin(); itLinker != mLinkers.end() && bIntegrity; ++itLinker)
		{
			bIntegrity = (itLinker->getStart().isValid() || itLinker->getEnd().isValid());
		}
		return bIntegrity;
	}
}

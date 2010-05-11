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
			if(1 == mLinkers.size())
			{
				const SecondaryStructure* pFront = mLinkers.front().start();
				const SecondaryStructure* pEnd = mLinkers.back().end();
				if(NULL != pFront && NULL != pEnd)
				{
					bClosed = pFront->isSame(*pEnd);
				}
			}
			else
			{
				bClosed = mLinkers.front().connects(mLinkers.back());
			}
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
			bool bStartConnected = (NULL != mLinkers.front().start());
			bool bEndConnected = (NULL != mLinkers.back().end());
			bComplete = (!bStartConnected && !bEndConnected);
		}
		return bComplete;
	}

	std::list<BaseInteraction> Loop::getBaseInteractions() const
	{
		std::list<BaseInteraction> interactions;
		std::vector<Linker>::const_iterator itLinker;
		std::vector<Linker>::const_iterator itNext;
		for(itLinker = mLinkers.begin(); itLinker != mLinkers.end(); ++itLinker)
		{
			itNext = itLinker;
			++ itNext;
			if(itNext == mLinkers.end())
			{
				itNext = mLinkers.begin();
			}

			std::list<BaseInteraction> linkerInteractions;
			linkerInteractions = itLinker->getBaseInteractions();
			interactions.insert(interactions.end(), linkerInteractions.begin(), linkerInteractions.end());
			std::list<BasePair> pairs = getLinkerPair(*itLinker, *itNext);
			if(0 < pairs.size())
			{
				if(pairs.front().fResId != interactions.back().rResId)
				{
					pairs.front().reverse();
				}
				interactions.push_back(pairs.front());
			}
		}
		return interactions;
	}


	std::list<BasePair> Loop::getLinkerPair(
		const Linker& aFirst,
		const Linker& aSecond) const
	{
		std::list<BasePair> pairs;
		if(aFirst.end() != NULL && aFirst.end() == aSecond.start())
		{
			// Both are pointing at the same stem
			const Stem* pStem = dynamic_cast<const Stem*>(aFirst.end());
			assert(NULL != pStem);
			// Note : Ignore faces of pair
			BasePair basePair(
				aFirst.residues().back().label(), aFirst.residues().back(),
				aSecond.residues().front().label(), aSecond.residues().front(),
				BasePair::face_vector());
			pairs.push_back(basePair);
		}
		return pairs;
	}

	std::set<mccore::ResId> Loop::getResIds() const
	{
		std::set<mccore::ResId> resids;
		std::vector<Linker>::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++ it)
		{
			resids.insert(it->residues().begin(), it->residues().end());
		}
		return resids;
	}

	bool Loop::checkIntegrity() const
	{
		bool bIntegrity = true;
		std::vector<Linker>::const_iterator itLinker;
		for(itLinker = mLinkers.begin(); itLinker != mLinkers.end() && bIntegrity; ++itLinker)
		{
			bIntegrity = (NULL != itLinker->start() || NULL != itLinker->end());
		}
		return bIntegrity;
	}

	std::set<LabeledResId> Loop::getSharedResIds() const
	{
		std::set<LabeledResId> resIds;
		std::vector<Linker>::const_iterator it;
		for(it = mLinkers.begin(); it != mLinkers.end(); ++ it)
		{
			std::set<LabeledResId> linkerResIds = it->getSharedResIds();
			resIds.insert(linkerResIds.begin(), linkerResIds.end());
		}
		return resIds;
	}
}

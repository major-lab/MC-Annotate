#include "Linker.h"

#include <cassert>

namespace annotate
{
	Linker::Linker()
	{
	}

	Linker::Linker(
		const std::vector<mccore::ResId>& aResidues,
		const StemConnection& aStart,
		const StemConnection& aEnd)
	{
		mResidues = aResidues;
		mStart = aStart;
		mEnd = aEnd;
	}

	Linker::~Linker()
	{
		clear();
	}

	void Linker::clear()
	{
		mResidues.clear();
	}

	bool Linker::isEmpty() const
	{
		bool bIsEmpty = true;

		if((mStart.isValid() && mEnd.isValid()) || 0 < mResidues.size())
		{
			bIsEmpty = false;
		}
		return bIsEmpty;
	}

	bool Linker::operator== (const Linker &other) const
	{
		bool bEqual = false;

		if((mStart == other.mStart && mEnd == other.mEnd) // same
			|| (mStart == other.mEnd && mEnd == other.mStart)) // reversed
		{
			// Same order
			bEqual = true;
		}
		return bEqual;
	}

	bool Linker::operator!= (const Linker &other) const
	{
		bool bEqual = operator==(other);

		return !bEqual;
	}

	bool Linker::operator< (const Linker &other) const
    {
    	bool bIsSmaller = false;

    	if(!isEmpty() && !other.isEmpty())
    	{
    		mccore::ResId thisId;
    		if(mStart.isValid())
	    	{
	    		thisId = mStart.getResidue();
	    	}
	    	else
	    	{
	    		thisId = mResidues.front();
	    	}

	    	mccore::ResId otherId;
	    	if(other.mStart.isValid())
	    	{
	    		otherId = other.mStart.getResidue();
	    	}
	    	else
	    	{
	    		otherId = other.mResidues.front();
	    	}
	    	bIsSmaller = thisId < otherId;
    	}

      return bIsSmaller;
    }

	bool Linker::isSame(const SecondaryStructure& aStruct) const
	{
		bool bSame = false;
		const Linker* pLinker = dynamic_cast<const Linker*>(&aStruct);
		if(NULL != pLinker && operator == (*pLinker))
		{
			// This is a linker and it has the same value as this one
			bSame = true;
		}
		return bSame;
	}

    bool Linker::isAdjacent(const SecondaryStructure& aStruct) const
	{
		bool bAdjacent = false;

		if(isSame(aStruct))
		{
			// Parameter is the current linker
			bAdjacent = true;
		}
		else if(mStart.isValid() && mStart.getStructure()->isSame(aStruct))
		{
			// Start of the linker connects to the provided structure
			bAdjacent = true;
		}
		else if(mEnd.isValid() && mEnd.getStructure()->isSame(aStruct))
		{
			bAdjacent = true;
		}
		return bAdjacent;
	}

    bool Linker::contains(const mccore::ResId& aResId) const
	{
		std::vector<mccore::ResId>::const_iterator it;
		it = std::find(mResidues.begin(), mResidues.end(), aResId);
		return it != mResidues.end();
	}

	void Linker::order()
	{
		if(0 < mResidues.size())
		{
			std::sort(mResidues.begin(), mResidues.end());
		}

		if(mStart.isValid() && mEnd.isValid())
		{
			mccore::ResId startResId = mStart.getResidue();
			mccore::ResId endResId = mEnd.getResidue();
			if(endResId < startResId)
			{
				std::swap(mStart, mEnd);
			}
		}else if(mStart.isValid() && 0 < mResidues.size())
		{
			mccore::ResId startResId = mStart.getResidue();
			if(startResId < mResidues.back())
			{
				std::swap(mStart, mEnd);
			}
		}else if(mEnd.isValid() && 0 < mResidues.size())
		{
			mccore::ResId endResId = mEnd.getResidue();
			if(endResId < mResidues.back())
			{
				std::swap(mStart, mEnd);
			}
		}
	}

	void Linker::reverse()
	{
		std::reverse(mResidues.begin(), mResidues.end());
		std::swap(mStart, mEnd);
	}

	bool Linker::connects(const Linker& aLinker) const
	{
		bool bConnects = false;
		if( mStart.connects(aLinker.mStart)
			|| mStart.connects(aLinker.mEnd)
			|| mEnd.connects(aLinker.mStart)
			|| mEnd.connects(aLinker.mEnd) )
		{
			bConnects = true;
		}
		return bConnects;
	}

	std::set<BaseInteraction> Linker::getBaseInteractions() const
	throw(mccore::FatalIntLibException)
	{
		std::set<BaseInteraction> interactions;
		std::vector<mccore::ResId>::const_iterator itRes;

		// Get the interactions between the residues of the linkers
		for(itRes = mResidues.begin(); itRes != mResidues.end(); ++itRes)
		{
			std::vector<mccore::ResId>::const_iterator itNextRes = itRes;
			itNextRes ++;
			if(itNextRes != mResidues.end())
			{
				// TODO : Interactions shouldn't use 0 labels, we need to get
				// the real adjacency interactions in the linkers
				BaseInteraction inter(0, *itRes, 0, *itNextRes);
				interactions.insert(inter);
			}
		}

		if(0 == mResidues.size())
		{
			appendInteractionBetweenConnections(interactions);
		}else
		{
			if(mStart.isValid())
			{
				interactions.insert(getBaseInteractionWithConnection(mStart));
			}
			if(mEnd.isValid())
			{
				interactions.insert(getBaseInteractionWithConnection(mEnd));
			}
		}

		return interactions;
	}

	void Linker::appendInteractionBetweenConnections(
			std::set<BaseInteraction>& aInteractionSet) const
			throw(mccore::FatalIntLibException)
	{
		if(mStart.isValid() && mEnd.isValid())
		{
			// Interaction is bewteen connections
			aInteractionSet.insert(getBaseInteractionBetweenConnections());
		}else
		{
			std::string strMsg("Linker::getBaseInteractions - Empty open linker");
			throw mccore::FatalIntLibException(strMsg, __FILE__, __LINE__);
		}
	}

	BaseInteraction Linker::getBaseInteractionWithConnection(
		const StemConnection& aConnection) const
		throw(mccore::FatalIntLibException)
	{
		assert(aConnection.isValid());
		assert(0 < mResidues.size());

		// TODO : Interactions shouldn't use 0 labels, we need to get
		// the real adjacency interactions in the linkers<
		res_info startResInfo = getResInfo(aConnection);
		res_info endResInfo;
		endResInfo.first = 0;

		if(startResInfo.second < std::min(mResidues.front(), mResidues.back()))
		{
			// Connection is 5'
			endResInfo.second = std::min(mResidues.front(), mResidues.back());
		}
		else
		{
			// Connection is 3'
			endResInfo.second = std::max(mResidues.front(), mResidues.back());
		}

		return BaseInteraction(
			startResInfo.first, startResInfo.second,
			endResInfo.first, endResInfo.second);
	}

	BaseInteraction Linker::getBaseInteractionBetweenConnections() const
		throw(mccore::FatalIntLibException)
	{
		// Assertions
		assert(mStart.isValid() && mEnd.isValid());

		res_info startResInfo = getResInfo(mStart);
		res_info endResInfo = getResInfo(mEnd);

		return BaseInteraction(
			startResInfo.first,
			startResInfo.second,
			endResInfo.first,
			endResInfo.second);
	}

	Linker::res_info Linker::getResInfo(const StemConnection& aConnection) const
		throw(mccore::NoSuchElementException)
	{
		std::pair<mccore::GraphModel::label, mccore::ResId> info;
		BasePair bp = aConnection.getPair();
		if(bp.fResId == aConnection.getResidue())
		{
			info.first = bp.first;
			info.second = bp.fResId;
		}else if(bp.rResId == aConnection.getResidue())
		{
			info.first = bp.second;
			info.second = bp.rResId;
		}else
		{
			std::string strMsg("Linker::getResInfo - Residue is not part of the connection");
			throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
		}
		return info;
	}

};

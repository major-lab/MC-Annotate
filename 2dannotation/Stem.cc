#include "Stem.h"

namespace annotate 
{
	void Stem::push_back(const BasePair& aBasePair)
	{
		// Insure that what we're trying to append is in order
		BasePair toAdd = orderPair(aBasePair);
		
		if(continues(toAdd))
		{
			if(1 == size())
			{
				const BasePair* pLast = &mBasePairs.back();
				meOrientation = pairOrientation(*pLast, toAdd);				
			}
			mBasePairs.push_back(toAdd);
			mResIds.insert(toAdd.fResId);
			mResIds.insert(toAdd.rResId);
		}
	}
	
	void Stem::clear()
	{
		meOrientation = eUNDEFINED;
		mBasePairs.clear();
		mResIds.clear();
	}
	
	bool Stem::operator ==(const Stem& other) const
	{
		bool bEqual = false;
		if(	mBasePairs.size() == other.mBasePairs.size() 
			&& meOrientation == other.meOrientation)
		{
			bEqual = std::equal(
				mBasePairs.begin(), 
				mBasePairs.end(), 
				other.mBasePairs.begin());
		}
		return bEqual;
	}
	
	bool Stem::contains(const mccore::Residue& aResidue) const
	{
		return contains(aResidue.getResId());
	}
	
	bool Stem::contains(const mccore::ResId& aResId) const
	{
		bool bContains = false;
		/*
		if(0 < size())
		{
			if(aResId == mBasePairs.front().fResId 
				|| aResId == mBasePairs.front().rResId 
				|| aResId == mBasePairs.back().fResId 
				|| aResId == mBasePairs.back().rResId)
			{
				bContains = true;
			}
			if(mBasePairs.front().fResId < aResId 
				&& aResId < mBasePairs.back().fResId)
			{
				bContains = true;
			}
			else if(aResId < std::max(mBasePairs.front().rResId, mBasePairs.back().rResId) 
				&& std::min(mBasePairs.front().rResId, mBasePairs.back().rResId) < aResId)
			{
				bContains = true;
			}
		}*/
		bContains = (mResIds.find(aResId) != mResIds.end());
		return bContains;
	}
	
	bool Stem::overlaps(const Stem& aStem) const
	{
		bool bOverlaps = false;
		if(0 < mBasePairs.size() && 0 < aStem.mBasePairs.size())
		{
			BasePair front = aStem.mBasePairs.front();
			BasePair back = aStem.mBasePairs.back();
			
			// Checks if we contains any of the bounds of the other
			bOverlaps = contains(front.fResId);
			bOverlaps |= contains(front.rResId);
			bOverlaps |= contains(back.fResId);
			bOverlaps |= contains(back.rResId);
			
			if(!bOverlaps)
			{
				// Check if the other contains this bounds
				bOverlaps = aStem.contains(mBasePairs.front().fResId);
				bOverlaps |= aStem.contains(mBasePairs.front().rResId);
				bOverlaps |= aStem.contains(mBasePairs.back().fResId);
				bOverlaps |= aStem.contains(mBasePairs.back().rResId);
			}
		}
		return bOverlaps;
	}
	
	bool Stem::continues(const BasePair& aBasePair) const
	{
		bool bContinues = true;
		
		if(0 < mBasePairs.size())
		{
			const BasePair* pLast = &mBasePairs.back();
			bContinues = pLast->fResId < aBasePair.fResId;
			bContinues = bContinues && pLast->areContiguous(aBasePair);
			if(1 < mBasePairs.size())
			{
				enOrientation eOrient;
				eOrient = pairOrientation(*pLast, aBasePair);
				bContinues = bContinues && (eOrient == meOrientation);
			}
		}
		return bContinues;
	}
	
	bool Stem::isAdjacent(const SecondaryStructure& aStruct) const
	{
		bool bAdjacent = false;
		
		const Stem* pStem = dynamic_cast<const Stem*>(&aStruct);
		if(NULL != pStem)
		{
			if(operator == (*pStem))
			{
				bAdjacent = true;
			}			
		}else
		{
			// Other structure is not a stem, ask it
			bAdjacent = aStruct.isAdjacent(*this);
		}
		return bAdjacent;
	}
	
	bool Stem::pseudoKnots( const Stem& aStem ) const
	{
		bool bPseudoKnots = false;
		if(0 < aStem.size() && 0 < size())
		{
			mccore::ResId aR1 = aStem.basePairs().front().fResId;
			mccore::ResId aR2 = aStem.basePairs().back().fResId;
			mccore::ResId aR3 = aStem.basePairs().front().rResId;
			mccore::ResId aR4 = aStem.basePairs().back().rResId;
			
			mccore::ResId tR1 = mBasePairs.front().fResId;
			mccore::ResId tR2 = mBasePairs.back().fResId;
			mccore::ResId tR3 = mBasePairs.front().rResId;
			mccore::ResId tR4 = mBasePairs.back().rResId;
			
			if(((tR2 < aR1 && aR2 < tR3) && (tR4 < aR3))
				|| ((aR2 < tR1 && tR2 < aR3) && (aR4 < tR3)))
			{
				bPseudoKnots = true;
			}
		}
		return bPseudoKnots;
	}
	
	std::vector<const Stem*> 
	Stem::pseudoKnots( std::vector<const Stem*>& aStems ) const
	{
		std::vector<const Stem*> oStems;
		std::vector<const Stem*>::const_iterator it = aStems.begin();
		for(;it != aStems.end(); ++it)
		{
			const Stem* pStem = *it;
			if(pseudoKnots(*pStem))
			{
				oStems.push_back(pStem);
			}
		}
		return oStems;
	}
	
	std::vector<const Stem*> 
	Stem::notPseudoKnotted( std::vector<const Stem*>& aStems ) const
	{
		std::vector<const Stem*> oStems;
		std::vector<const Stem*>::const_iterator it = aStems.begin();
		for(;it != aStems.end(); ++it)
		{
			const Stem* pStem = *it;
			if(!pseudoKnots(*pStem))
			{
				oStems.push_back(pStem);
			}
		}
		return oStems;
	}
	
	BasePair Stem::orderPair(const BasePair& aBasePair) const
	{
		// Insure that what we're trying to append is in order
		BasePair orderedPair = aBasePair;
		if(orderedPair.rResId < orderedPair.fResId)
		{
			orderedPair.reverse();
		}
		return orderedPair;
	}

	Stem::enOrientation Stem::pairOrientation(
		const BasePair& aBasePair1, 
		const BasePair& aBasePair2)
	{
		enOrientation eOrientation = eUNDEFINED;
		int iR1 = aBasePair1.rResId.getResNo() - aBasePair2.rResId.getResNo();
		int iR2 = aBasePair1.fResId.getResNo() - aBasePair2.fResId.getResNo();
		
		if((iR1 < 0 && iR2 < 0) || (iR1 > 0 && iR2 > 0))
		{
			eOrientation = ePARALLEL;
		}
		else if((iR1 < 0 && iR2 > 0) || (iR1 > 0 && iR2 < 0))
		{
			eOrientation = eANTIPARALLEL;
		}
		return eOrientation;
	}
	
	mccore::ResId Stem::getResidue(const enConnection& ePoint) const 
		throw(mccore::NoSuchElementException)
	{
		mccore::ResId residueId;
		switch(ePoint)
		{
			case eFIRST_STRAND_FRONT_PAIR:
				residueId = mBasePairs.front().fResId;
				break;
			case eSECOND_STRAND_FRONT_PAIR:
				residueId = mBasePairs.front().rResId;
				break;
			case eFIRST_STRAND_BACK_PAIR:
				residueId = mBasePairs.back().fResId;
				break;	
			case eSECOND_STRAND_BACK_PAIR:
				residueId = mBasePairs.back().rResId;
				break;
			default:
				std::string strMsg("ePoint must be defined");
				throw mccore::NoSuchElementException(strMsg, __FILE__, __LINE__);
		}
		return residueId;
	}
	
	Stem::enConnection Stem::getConnection(const mccore::ResId& aId) const
	{
		enConnection eConnect = eUNDEFINED_CONNECTION;
		
		if(0 < mBasePairs.size())
		{
			if(mBasePairs.front().fResId == aId)
			{
				eConnect = eFIRST_STRAND_FRONT_PAIR;
			}
			else if(mBasePairs.front().rResId == aId)
			{
				eConnect = eSECOND_STRAND_FRONT_PAIR;
			}
			else if(mBasePairs.back().fResId == aId)
			{
				eConnect = eFIRST_STRAND_BACK_PAIR;
			}
			else if(mBasePairs.back().rResId == aId)
			{
				eConnect = eSECOND_STRAND_BACK_PAIR;
			}
		}
		return eConnect;
	}
	
	int StemConnection::getDirection() const
	{
		int iDirection = 0;
		if(isValid())
		{
			switch(meConnection)
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
				iDirection = -1;
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
				iDirection = 1;
				break;
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				if(Stem::eANTIPARALLEL == mpStem->getOrientation())
				{
					iDirection = 1;
				}
				else
				{
					iDirection = -1;
				}
				break;
			case Stem::eSECOND_STRAND_BACK_PAIR:
				if(Stem::eANTIPARALLEL == mpStem->getOrientation())
				{
					iDirection = -1;
				}
				else
				{
					iDirection = 1;
				}
				break;
			default:
				// TODO : Exception
				iDirection = 0;
			}
		}
		return iDirection;
	}
	
	mccore::ResId StemConnection::nextId() const
	{
		mccore::ResId id;
		if(NULL != mpStem && meConnection != Stem::eUNDEFINED_CONNECTION)
		{
			switch(meConnection)
			{
			case Stem::eFIRST_STRAND_FRONT_PAIR:
				id = mpStem->basePairs().front().rResId;
				break;
			case Stem::eSECOND_STRAND_FRONT_PAIR:
				id = mpStem->basePairs().front().fResId;
				break;
			case Stem::eFIRST_STRAND_BACK_PAIR:
				id = mpStem->basePairs().back().rResId;
				break;
			case Stem::eSECOND_STRAND_BACK_PAIR:
				id = mpStem->basePairs().back().fResId;
				break;
			default:
				// TODO : This should be an exception
				break;
			}
		}
		else
		{
			// TODO : Throw an exception
		}
		
		return id;
	}
};


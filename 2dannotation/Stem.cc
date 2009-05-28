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
		}
	}
	
	void Stem::clear()
	{
		meOrientation = eUNDEFINED;
		mBasePairs.clear();
	}
		
	bool Stem::continues(const BasePair& aBasePair) const
	{
		bool bContinues = true;
		
		if(0 < mBasePairs.size())
		{
			const BasePair* pLast = &mBasePairs.back();
			bContinues = pLast->fResId < aBasePair.fResId;
			bContinues &= pLast->areContiguous(aBasePair);
			if(1 < mBasePairs.size())
			{
				enOrientation eOrient;
				eOrient = pairOrientation(*pLast, aBasePair);
				bContinues &= (eOrient == meOrientation);
			}
		}
		return bContinues;
	}
	
	bool Stem::pseudoKnots( const Stem& aStem ) const
	{
		bool bPseudoKnots = false;
		if(0 < aStem.size() && 0 < size())
		{
			ResId aR1 = aStem.basePairs().front().fResId;
			ResId aR2 = aStem.basePairs().back().fResId;
			ResId aR3 = aStem.basePairs().front().rResId;
			ResId aR4 = aStem.basePairs().back().rResId;
			
			ResId tR1 = mBasePairs.front().fResId;
			ResId tR2 = mBasePairs.back().fResId;
			ResId tR3 = mBasePairs.front().rResId;
			ResId tR4 = mBasePairs.back().rResId;
			
			if(((tR2 < aR1 && aR2 < tR3) && (tR4 < aR3))
				|| ((aR2 < tR1 && tR2 < aR3) && (aR4 < tR3)))
			{
				bPseudoKnots = true;
			}
		}
		return bPseudoKnots;
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
};


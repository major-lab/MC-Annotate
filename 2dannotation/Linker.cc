#include "Linker.h"
#include "BaseLink.h"

#include <cassert>

namespace annotate
{
	Linker::Linker(
		const std::vector<LabeledResId>& aResidues,
		const SecondaryStructure* apStartStruct,
		const SecondaryStructure* apEndStruct)
	{
		mResidues = aResidues;
		mpStartStruct = apStartStruct;
		mpEndStruct = apEndStruct;
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
		return (2 == mResidues.size());
	}

	bool Linker::operator== (const Linker &other) const
	{
		bool bEqual = false;

		if(mResidues.size() == other.mResidues.size())
		{
			if(0 < mResidues.size())
			{
				bEqual = (mResidues.front() == other.mResidues.front()
					&& mResidues.back() == other.mResidues.back());
			}else
			{
				bEqual = true;
			}
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
	    	bIsSmaller = mResidues.front() < other.mResidues.front();
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
		else if(NULL != mpStartStruct && mpStartStruct->isSame(aStruct))
		{
			// Start of the linker connects to the provided structure
			bAdjacent = true;
		}
		else if(NULL != mpEndStruct && mpEndStruct->isSame(aStruct))
		{
			bAdjacent = true;
		}
		return bAdjacent;
	}

    bool Linker::contains(const LabeledResId& aResId) const
	{
		std::vector<LabeledResId>::const_iterator it;
		it = std::find(mResidues.begin(), mResidues.end(), aResId);
		return it != mResidues.end();
	}

	void Linker::order()
	{
		assert(2 <= mResidues.size());

		if(mResidues.back() < mResidues.front())
		{
			reverse();
		}
	}

	void Linker::reverse()
	{
		std::reverse(mResidues.begin(), mResidues.end());
		std::swap(mpStartStruct, mpEndStruct);
	}

	bool Linker::connects(const Linker& aLinker) const
	{
		bool bConnects = false;
		if( mResidues.front() == aLinker.mResidues.front()
			|| mResidues.front() == aLinker.mResidues.back()
			|| mResidues.back() == aLinker.mResidues.front()
			|| mResidues.back() == aLinker.mResidues.back())
		{
			bConnects = true;
		}else if(NULL != mpStartStruct
					&& (mpStartStruct == aLinker.mpStartStruct
						|| mpStartStruct == aLinker.mpEndStruct))
		{
			bConnects = true;
		}else if(NULL != mpEndStruct
				&& (mpEndStruct == aLinker.mpStartStruct
					|| mpEndStruct == aLinker.mpEndStruct))
		{
			bConnects = true;
		}

		return bConnects;
	}

	std::set<BaseInteraction> Linker::getBaseInteractions() const
	throw(mccore::FatalIntLibException)
	{
		std::set<BaseInteraction> interactions;
		std::vector<LabeledResId>::const_iterator itRes;

		// Get the interactions between the residues of the linkers
		for(itRes = mResidues.begin(); itRes != mResidues.end(); ++itRes)
		{
			std::vector<LabeledResId>::const_iterator itNextRes = itRes;
			itNextRes ++;
			if(itNextRes != mResidues.end())
			{
				BaseLink inter(
					itRes->label(), *itRes,
					itNextRes->label(), *itNextRes);
				interactions.insert(inter);
			}
		}

		return interactions;
	}

	std::set<LabeledResId> Linker::getSharedResIds() const
	{
		std::set<LabeledResId> resIds;
		if(NULL != mpStartStruct)
		{
			resIds.insert(mResidues.front());
		}
		if(NULL != mpEndStruct)
		{
			resIds.insert(mResidues.back());
		}
		return resIds;
	}

};

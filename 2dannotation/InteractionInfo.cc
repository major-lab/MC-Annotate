#include "InteractionInfo.h"

#include <cassert>

namespace annotate {

InteractionInfo::InteractionInfo(
		const std::string& astrFile,
		unsigned int auiModel,
		const mccore::ResId& aResId1,
		const mccore::ResId& aResId2,
		const enFace& aFace1,
		const enFace& aFace2,
		const annotate::BasePair::enOrientation& aeOrientation)
: mModelInfo(astrFile, auiModel)
{
	// TODO : Insure that the residues are in 5' -> 3' order
	assert(aResId1 < aResId2);
	mResId1 = aResId1;
	mResId2 = aResId2;
	meFace1 = aFace1;
	meFace2 = aFace2;
	meOrientation = aeOrientation;
}

InteractionInfo::InteractionInfo(
		const std::string& astrFile,
		unsigned int auiModel,
		const std::string& astrRes1,
		const std::string& astrRes2,
		const enFace& aFace1,
		const enFace& aFace2,
		const annotate::BasePair::enOrientation& aeOrientation)
: mModelInfo(astrFile, auiModel)
{
	// TODO : Insure that the residues are in 5' -> 3' order
	mResId1 = mccore::ResId(astrRes1.c_str());
	mResId2 = mccore::ResId(astrRes2.c_str());
	assert(mResId1 < mResId2);

	meFace1 = aFace1;
	meFace2 = aFace2;
	meOrientation = aeOrientation;
}

bool InteractionInfo::operator <(const InteractionInfo& aRight) const
{
	bool bLess = false;

	if(mModelInfo < aRight.mModelInfo)
	{
		bLess = true;
	}
	else if(mModelInfo == aRight.mModelInfo)
	{
		if(mResId1 < aRight.mResId1)
		{
			bLess = true;
		} else if(mResId1 == aRight.mResId1)
		{
			if(mResId2 < aRight.mResId1)
			{
				bLess = true;
			}else if(mResId2 < aRight.mResId1)
			{
				if(meFace1 < aRight.meFace1)
				{
					bLess = true;
				}else if(meFace1 == aRight.meFace1)
				{
					if(meFace2 < aRight.meFace2)
					{
						bLess = true;
					}else if(meFace2 == aRight.meFace2)
					{
						bLess = (meOrientation < aRight.meOrientation);
					}
				}
			}
		}
	}
	return bLess;
}

bool InteractionInfo::operator ==(const InteractionInfo& aRight) const
{
	bool bEqual = (		(mModelInfo == aRight.mModelInfo)
		&&	(mResId1 == aRight.mResId1)
		&&  (mResId2 == aRight.mResId2)
		&&  (meFace1 == aRight.meFace1)
		&&  (meFace2 == aRight.meFace2)
		&&  (meOrientation == aRight.meOrientation));
	return bEqual;
}

std::set<char> InteractionInfo::getChains() const
{
	std::set<char> chains;
	chains.insert(mResId1.getChainId());
	chains.insert(mResId2.getChainId());
	return chains;
}

void InteractionInfo::setChainAndOffset(char acChain, int aiOffset)
{
	mResId1.setChainId(acChain);
	mResId1.setResNo(mResId1.getResNo() + aiOffset);
	mResId2.setChainId(acChain);
	mResId2.setResNo(mResId2.getResNo() + aiOffset);
}

}; // namespace annotate

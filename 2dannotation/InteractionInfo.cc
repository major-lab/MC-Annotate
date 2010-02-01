#include "InteractionInfo.h"

namespace annotate {

InteractionInfo::InteractionInfo(
	const std::string& astrFile,
	unsigned int auiModel,
	const std::string& astrRes1,
	const std::string& astrRes2)
: mModelInfo(astrFile, auiModel)
{
	// TODO : Insure that the residues are in 5' -> 3' order
	mstrRes1 = astrRes1;
	mstrRes2 = astrRes2;
}

InteractionInfo& InteractionInfo::operator=(const InteractionInfo& aRight)
{
	return *this;
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
		if(mstrRes1 < aRight.mstrRes1)
		{
			bLess = true;
		} else if(mstrRes1 == aRight.mstrRes1)
		{
			if(mstrRes2 < aRight.mstrRes2)
			{
				bLess = true;
			}
		}
	}
	return bLess;
}

}; // namespace annotate

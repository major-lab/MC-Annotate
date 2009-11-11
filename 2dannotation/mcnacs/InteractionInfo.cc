#include "InteractionInfo.h"

InteractionInfo::InteractionInfo(
	const std::string& astrFile,
	unsigned int auiModel,
	const std::string& astrRes1,
	const std::string& astrRes2)
: mModelInfo(astrFile, auiModel)
{
	mstrRes1 = astrRes1;
	mstrRes2 = astrRes2;
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

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
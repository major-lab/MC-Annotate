#include "ModelInfo.h"

namespace annotate {

ModelInfo::ModelInfo(const std::string& astrFile, unsigned int auiModel)
{
	mstrPDBFile = astrFile;
	muiModel = auiModel;
}

bool ModelInfo::operator <(const ModelInfo& aRight) const
{
	return ((mstrPDBFile < aRight.mstrPDBFile)
		|| (mstrPDBFile == aRight.mstrPDBFile && (muiModel < aRight.muiModel)));
}

bool ModelInfo::operator ==(const ModelInfo& aRight) const
{
	return ((muiModel == aRight.muiModel) && (mstrPDBFile == aRight.mstrPDBFile));
}

bool ModelInfo::operator !=(const ModelInfo& aRight) const
{
	return !operator ==(aRight);
}

}; // namespace annotate

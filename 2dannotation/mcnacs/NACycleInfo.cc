#include "NACycleInfo.h"

#include <cassert>

NACycleInfo::NACycleInfo(
	const std::string& aFile, 
	unsigned int auiModel, 
	const CycleInfo::residue_profile& aResidues) 
: CycleInfo(aFile, auiModel, aResidues)
{
	mConnections.resize(aResidues.size());
}

NACycleInfo::NACycleInfo(const CycleInfo& aCycle)
: CycleInfo(aCycle)
{
	mConnections.resize(mResidues.size());	
}

std::set<CycleInfo>& NACycleInfo::getStrandConnections(unsigned int auiStrand)
{
	assert(auiStrand < mConnections.size());
	return mConnections[auiStrand];
}
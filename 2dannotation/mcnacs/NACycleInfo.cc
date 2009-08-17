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

bool NACycleInfo::operator <(const NACycleInfo& aRight) const
{
	bool bLess = false;
	
	if(CycleInfo::operator <(aRight))
	{
		bLess = true;
	}else
	{
		const CycleInfo* pRight = dynamic_cast<const CycleInfo*>(&aRight);
		assert(NULL != pRight);
		if(pRight->operator <(*this))
		{
			bLess = false;
		}
		else
		{
			// Equal at cycle info level, must check the connected cycles
			for(unsigned int i = 0; i < mConnections.size(); ++ i)
			{
				int iCompare = compareConnectionStrand(
					mConnections[i], 
					aRight.mConnections[i]);
				if(iCompare < 0)
				{
					bLess = true;
					break;
				}else if(iCompare > 0)
				{
					bLess = false;
					break;
				}
			}
		}
	}
	return bLess;
}

int NACycleInfo::compareConnectionStrand(
	const std::set<CycleInfo>& aConnection1, 
	const std::set<CycleInfo>& aConnection2) const
{
	int iCompare = 0;
	
	if(aConnection1.size() < aConnection2.size())
	{
		iCompare = -1;
	}else if(aConnection2.size() < aConnection1.size())
	{
		iCompare = 1;
	}else
	{
		std::set<CycleInfo>::const_iterator it1 = aConnection1.begin();
		std::set<CycleInfo>::const_iterator it2 = aConnection2.begin();
		while(it1 != aConnection1.end() && it2 != aConnection2.end())
		{
			if(*it1 < *it2)
			{
				iCompare = -1;
				break;
			}else if(*it2 < *it1)
			{
				iCompare = 1;
				break;
			}
			++ it1;
			++ it2;
		}
	}
	return iCompare;
}
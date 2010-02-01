#include "NACycleInfo.h"

#include "CycleProfile.h"

#include <cassert>
#include <sstream>

NACycleInfo::NACycleInfo(
	const std::string& aFile,
	unsigned int auiModel,
	const annotate::CycleProfile& aFileProfile,
	const annotate::CycleProfile& aProfile,
	const annotate::CycleInfo::residue_profile& aResIds,
	const std::vector<std::string>& aResidues)
: annotate::CycleInfo(aFile, auiModel, aFileProfile, aProfile, aResIds, aResidues)
{
	mConnections.resize(aResIds.size());
}

NACycleInfo::NACycleInfo(const annotate::CycleInfo& aCycle)
: annotate::CycleInfo(aCycle)
{
	mConnections.resize(mResIds.size());
}

std::set<annotate::CycleInfo>& NACycleInfo::getStrandConnections(
	unsigned int auiStrand)
{
	assert(auiStrand < mConnections.size());
	return mConnections[auiStrand];
}

bool NACycleInfo::operator <(const NACycleInfo& aRight) const
{
	bool bLess = false;

	if(annotate::CycleInfo::operator <(aRight))
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

std::string NACycleInfo::equationString() const
{
	std::ostringstream oss;
	std::string strProfile = this->getProfile().toString();
	oss << strProfile << " = ";
	std::vector< std::set<CycleInfo> >::const_iterator it;
	for(it = mConnections.begin(); it != mConnections.end(); ++ it)
	{
		if(it != mConnections.begin())
		{
			oss << " + ";
		}
		if(1 < it->size())
		{
			oss << "(";
		}
		std::set<CycleInfo>::const_iterator itCycle;
		for(itCycle = it->begin(); itCycle != it->end(); ++ itCycle)
		{
			if(itCycle != it->begin())
			{
				oss << " + ";
			}
			oss << itCycle->getProfile().toString();
		}
		if(1 < it->size())
		{
			oss << ")";
		}
	}
	return oss.str();
}

std::string NACycleInfo::modelInfoString() const
{
	std::ostringstream oss;
	oss << this->getModelInfo().getPDBFile();
	oss << "[" << this->getModelInfo().getModel() << "]";
	return oss.str();
}

std::string NACycleInfo::toString() const
{
	std::ostringstream oss;
	oss << "index:" << std::endl;
	oss << "\t" << equationString() << std::endl;
	oss << "info:" << std::endl;
	oss << "\t" << modelInfoString() << std::endl;
	oss << "\t" << CycleInfo::toString() << std::endl;
	oss << "\t=" << std::endl;
	std::vector<std::set<CycleInfo> >::const_iterator it;
	for(it = mConnections.begin(); it != mConnections.end(); ++ it)
	{
		if(it != mConnections.begin())
		{
			oss << "\t+" << std::endl;
		}
		std::set<CycleInfo>::const_iterator itCycle;
		for(itCycle = it->begin(); itCycle != it->end(); ++itCycle)
		{
			oss << "\t" << itCycle->toString() << std::endl;
		}
	}
	return oss.str();
}

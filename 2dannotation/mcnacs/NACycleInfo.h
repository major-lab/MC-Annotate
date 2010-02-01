#ifndef _NACycleInfo_H_
#define _NACycleInfo_H_

#include "CycleInfo.h"
#include <set>

namespace annotate
{
	class CycleProfile;
};

class NACycleInfo : public annotate::CycleInfo
{
public:
	// LIFECYLE ----------------------------------------------------------------
	NACycleInfo(
		const std::string& aFile,
		unsigned int auiModel,
		const annotate::CycleProfile& aFileProfile,
		const annotate::CycleProfile& aProfile,
		const annotate::CycleInfo::residue_profile& aResIds,
		const std::vector<std::string>& aResidues);
	NACycleInfo(const annotate::CycleInfo& aCycle);
	~NACycleInfo(){}

	// ACCESS ------------------------------------------------------------------
	std::set<annotate::CycleInfo>& getStrandConnections(unsigned int auiStrand);
	const std::vector<std::set<annotate::CycleInfo> >& getConnections() const
	{return mConnections;}

	// OPERATOR ----------------------------------------------------------------
	bool operator <(const NACycleInfo& aRight) const;

	// METHODS -----------------------------------------------------------------
	std::string toString() const;
	std::string equationString() const;

private:
	std::vector<std::set<annotate::CycleInfo> > mConnections;

	int compareConnectionStrand(
		const std::set<annotate::CycleInfo>& aConnection1,
		const std::set<annotate::CycleInfo>& aConnection2) const;
	std::string modelInfoString() const;
};

#endif /*_NACycleInfo_H_*/

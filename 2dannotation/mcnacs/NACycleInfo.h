#ifndef _NACycleInfo_H_
#define _NACycleInfo_H_

#include "CycleInfo.h"

#include <set>

class NACycleInfo : public CycleInfo
{
public:
	// LIFECYLE ----------------------------------------------------------------
	NACycleInfo(
		const std::string& aFile, 
		unsigned int auiModel, 
		const CycleInfo::residue_profile& aResidues);
	NACycleInfo(const CycleInfo& aCycle);
	~NACycleInfo(){}
	
	// ACCESS ------------------------------------------------------------------
	std::set<CycleInfo>& getStrandConnections(unsigned int auiStrand);
	const std::vector< std::set<CycleInfo> >& getConnections() const 
	{return mConnections;}
	
	// OPERATOR ----------------------------------------------------------------
	bool operator <(const NACycleInfo& aRight) const;
	
private:
	std::vector< std::set<CycleInfo> > mConnections;
	
	int compareConnectionStrand(
		const std::set<CycleInfo>& aConnection1, 
		const std::set<CycleInfo>& aConnection2) const;
};

#endif /*_NACycleInfo_H_*/
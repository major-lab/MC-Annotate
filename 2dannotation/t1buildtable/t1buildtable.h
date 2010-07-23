#ifndef _t1buildtable_h_
#define _t1buildtable_h_

#include "CycleInfo.h"
#include "RNAGroupFile.h"
#include "RNAGroupModel.h"

#include <map>
#include <set>
#include <string>
#include <vector>

class T1BuildTable
{
public:
	// LIFECYCLE
	T1BuildTable(int argc, char * argv []);

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;
private:
	typedef std::set<annotate::CycleInfo> cycle_set;
	void readOptions(int argc, char* argv[]);

	void readCyclesCount();
	void readPDBGroups();
	void readCyclesPairs();
	void computeNormalizedFrequencies();
	unsigned int profileId(const std::string& astrProfile) const;
	void outputFrequencies() const;

	std::string mstrCyclesCountFile;
	std::string mstrPDBGroupsFile;
	std::string mstrCyclesPairsFile;

	std::map<std::string, unsigned int> mProfileMap;

	std::map<unsigned int, annotate::RNAGroupModel> mGroupModelMap;

	std::vector<unsigned int> mCycleCount;
	annotate::RNAGroupFile mGroupFile;
	std::vector<std::vector<unsigned int> > mCyclesPairsCount;
	std::vector<std::vector<float> > mInteractionFrequencies;

	annotate::CycleInfo getCycleInfo(
		const std::string& astrPDB, unsigned int uiModel,
		const std::string& astrProfile,
		const std::string& astrResIds,
		const std::string& astrSequence) const;

	annotate::CycleInfo::residue_strand getResIds(
		const std::string& aResidues) const;
};

#endif /*_t1buildtable_h_*/

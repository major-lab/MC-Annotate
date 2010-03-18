#ifndef _mc3dicp_h_
#define _mc3dicp_h_

#include "CycleInfo.h"
#include <string>

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
	void readCyclesPairs();
	void computeNormalizedFrequencies();
	unsigned int profileId(const std::string& astrProfile) const;
	void outputFrequencies() const;

	std::string mstrCyclesCountFile;
	std::string mstrCyclesPairsFile;

	std::map<std::string, unsigned int> mProfileMap;
	std::vector<unsigned int> mCycleCount;
	std::vector<std::vector<unsigned int> > mCyclesPairsCount;
	std::vector<std::vector<float> > mInteractionFrequencies;
};

#endif /*_mc3dicp_h_*/

#ifndef _db2ncm_h_
#define _db2ncm_h_

#include "CycleInfo.h"

#include <map>
#include <set>
#include <string>

class DotBracketToNCM
{
public:
	// LIFECYCLE
	DotBracketToNCM(int argc, char * argv []);

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;
private:
	typedef std::pair<unsigned int, unsigned int> index_pair;
	void readOptions(int argc, char* argv[]);

	double mfProbability;

	std::string mstrIdentifier;
	std::string mstrSequence;
	std::string mstrDotBracket;

	std::vector<index_pair> mPairs;
	std::set<annotate::CycleInfo> mNCMs;

	std::vector<index_pair> identifyPairs() const;
	std::set<annotate::CycleInfo> identifyNCMs(
		unsigned int auiDBId,
		const std::vector<index_pair>& aPairs) const;
	// Double stranded
	annotate::CycleInfo getCycleInfo(
		unsigned int auiDBId,
		unsigned int auiCycleId,
		const std::vector<mccore::ResId>& aStrand1,
		const std::vector<mccore::ResId>& aStrand2) const;
	// Loop
	annotate::CycleInfo getCycleInfo(
		unsigned int auiDBId,
		unsigned int auiCycleId,
		const std::vector<mccore::ResId>& aStrand1) const;

	std::string cycleName(
		unsigned int auiDBId,
		unsigned int auiCycleId) const;
};

#endif /*_db2ncm_h_*/

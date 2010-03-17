#ifndef _mc3dicp_h_
#define _mc3dicp_h_

#include "CycleInfo.h"
#include "InteractionInfo.h"

#include <set>
#include <string>

class MC3DInteractingCyclePairs
{
public:
	// LIFECYCLE
	MC3DInteractingCyclePairs(int argc, char * argv []);

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;
private:
	typedef std::set<annotate::CycleInfo> cycle_set;
	typedef std::set<annotate::InteractionInfo> interaction_set;
	void readOptions(int argc, char* argv[]);

	std::string mstrCyclesFile;
	std::string mstrInteractionsFile;
	std::map<annotate::ModelInfo, cycle_set > mCycles;
	std::map<annotate::ModelInfo, interaction_set > mInteractions;
	std::set<annotate::ModelInfo> mModels;

	void identifyPairs();
	void identifyModelPairs(const annotate::ModelInfo& aModel);
	void outputPairs(
		std::set<annotate::CycleInfo>& aLeft,
		std::set<annotate::CycleInfo>& aRight);
};

#endif /*_mc3dicp_h_*/

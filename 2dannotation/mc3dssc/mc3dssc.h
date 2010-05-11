#ifndef _mc3dssc_h_
#define _mc3dssc_h_

#include "CycleInfo.h"
#include "InteractionInfo.h"

#include <set>
#include <string>

class MC3DSecondaryStructureCycles
{
public:
	// LIFECYCLE
	MC3DSecondaryStructureCycles(int argc, char * argv []);

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

	void identifyCycles();
	void identifyModelCycles(const annotate::ModelInfo& aModel);
	void outputCycles(std::set<annotate::CycleInfo>& aCycles) const;
	bool singleChain(const annotate::CycleInfo& aCycle) const;
};

#endif /*_mc3dssc_h_*/

#ifndef _mc3di_h_
#define _mc3di_h_

#include "CycleInfo.h"
#include "InteractionInfo.h"

#include <set>
#include <string>

class MC3DInteractions
{
public:
	// LIFECYCLE
	MC3DInteractions(int argc, char * argv []);

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
	void outputInteractions(
		const annotate::InteractionInfo& aInteraction,
		const std::set<annotate::CycleInfo>& aLeft,
		const std::set<annotate::CycleInfo>& aRight);
	std::string faceString(
		const annotate::InteractionInfo::enFace& aeFace) const;
	std::string orientationString(
		const annotate::BasePair::enOrientation& aeOrientation) const;
};

#endif /*_mc3dicp_h_*/

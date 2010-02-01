/*
 * NonAdjacentCycleStatsApp.h
 *
 *  Created on: Dec 15, 2009
 *      Author: blanchmf
 */

#ifndef NONADJACENTCYCLESTATSAPP_H_
#define NONADJACENTCYCLESTATSAPP_H_

#include "CycleInfo.h"
#include "CycleStatsEntry.h"
#include "InteractionInfo.h"
#include "InteractionTable.h"
#include "GetModelRangeFunctor.h"
#include "NACycleInfo.h"
#include "NAInteractionInfo.h"
#include <set>
#include <string>

class NonAdjacentCycleStatsApp
{
public:
	/**
	 * Constructor.
	 */
	NonAdjacentCycleStatsApp(int argc, char * argv []);

	// Console specific methods
	void version() const;
	void usage() const;
	void help() const;

	// Temporary
	bool removeComposites() const {return mbRemoveComposites;}
	bool splitInteractions() const {return mbSplitInteractions;}
	const std::set<annotate::InteractionInfo>& interactions() const {return mInteractions;}
	const std::set<annotate::CycleInfo>& cycles() const {return mCycles;}
	const std::set<annotate::CycleInfo>& secondaryStructureCycles() const {return mSecondaryStructureCycles;}
	const std::set<NACycleInfo>& connectedCycles() const {return mConnectedCycles;}
	void connectedCycles(const std::set<NACycleInfo>& aCycles) {mConnectedCycles = aCycles;}
	std::string toString() const;

	// METHODS -----------------------------------------------------------------
	std::set<annotate::CycleInfo> getNonAdjacentCycles() const;
	std::set<annotate::CycleInfo> getDoubleInteractionCycles() const;
	void computeCyclesStats();
	/*
	 * @brief Remove all connections that doesn't form an NCM on one side
	 */
	void filterOutOpenConnections();
	void filterOutPartialCoverage();
	void filterOutEnclosingCycles();
	void filterOutMultibranchCycles();
	void filterOutLooseCycles();
	void splitAdjacency();
	void compileStatistics();

private:
	bool mbSplitInteractions;
	bool mbRemoveComposites;
	std::string mstrCyclesFile;
	std::string mstrInteractionsFile;
	std::string mstrOutputDirectory;
	std::set<annotate::InteractionInfo> mInteractions;
	std::set<annotate::CycleInfo> mCycles;
	std::set<annotate::CycleInfo> mNACycles;
	std::set<annotate::CycleInfo> mSecondaryStructureCycles;
	std::set<NACycleInfo> mConnectedCycles; // May need to be removed
	std::map<std::string, unsigned int> mStatistics;
	std::set<CycleStatsEntry, lessCycleStatsEntry> mCyclesStats;
	InteractionTable mInteractionTable;
	unsigned int muiOpenConnections;
	unsigned int muiPartialConnections;
	unsigned int muiEnclosingCycles; // Number of cycles with enclosed cycles
	unsigned int muiMultibranchCycles; // Number of multibranch cycles
	unsigned int muiLooseCycles; // Number of loose connection cycles

	void readOptions(int argc, char* argv[]);
	void readInteractionsFile();
	void readCyclesFile();
	std::vector<std::vector<std::string> > getStrandResidues(
		const std::string& aResidues,
		const annotate::CycleProfile& aProfile) const;
	std::list<std::string> getResidues(const std::string& aResidues) const;

	std::set<annotate::CycleInfo> getNonAdjacentCycleFromModel(
		const GetModelRangeFunctor<annotate::InteractionInfo>::const_range& interactionRange,
		const GetModelRangeFunctor<annotate::CycleInfo>::const_range& cycleRange) const;
	std::list<annotate::ModelInfo> getModels() const;
	std::set<annotate::CycleInfo> getCyclesWithInteraction(
		const GetModelRangeFunctor<annotate::CycleInfo>::const_range aRange,
		const std::set<annotate::Interaction>& aInteractions) const;
	std::set<NACycleInfo> getConnectedCycles(
		const std::set<annotate::CycleInfo>& aNACycles) const;

	std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > >
	getDoubleInteractionCyclesPairsFromModel(
		const std::set<NAInteractionInfo>& aInteractions,
		const annotate::ModelInfo& aModel) const;
	std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > >
	getDoubleInteractionCyclesPairs(
		const std::set<annotate::InteractionInfo>& aInteractions,
		const std::set<annotate::CycleInfo>& aCycles) const;
	std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > >
	getDoubleInteractionCyclesPairsFromInteractionRange(
		const GetModelRangeFunctor<NAInteractionInfo>::const_range& aInterRange) const;
	std::set<CycleStatsEntry, lessCycleStatsEntry> computeCycleStatsInstances() const;
	void addInteractionBetween(std::set<CycleStatsEntry, lessCycleStatsEntry>& aCycleStats,
		const annotate::CycleInfo& aCycle1,
		const annotate::CycleInfo& aCycle2);
	void computeCycleStatsInteractions(
		std::set<CycleStatsEntry, lessCycleStatsEntry>& aCycleStats,
		const std::set<annotate::InteractionInfo>& aInteractions,
		const std::set<annotate::CycleInfo>& aCycles);
	std::set<std::pair< annotate::CycleInfo, annotate::CycleInfo > >
	getInteractingCyclesPairsFromInteraction(const NAInteractionInfo& aInteraction) const;
	std::set<NAInteractionInfo> getNAInteractions(
		const std::set<annotate::InteractionInfo>& aInteractions,
		const std::set<annotate::CycleInfo>& aCycles) const;
	std::set<NAInteractionInfo> getNAInteractionsFromModel(
		const GetModelRangeFunctor<annotate::InteractionInfo>::const_range& interactionRange,
		const GetModelRangeFunctor<annotate::CycleInfo>::const_range& cycleRange) const;
	std::string statisticsToString() const;
	std::string cycleStatsToString() const;
	std::string interactionsToString() const;
	std::string interactionStats(
		const InteractionTable::interacting_set& aInteracting) const;
	void outputInteractingCyclesFiles(
		const std::string& astrProfile1,
		const std::string& astrProfile2,
		const InteractionTable::interacting_set& aInteracting) const;
};

#endif /* NONADJACENTCYCLESTATSAPP_H_ */

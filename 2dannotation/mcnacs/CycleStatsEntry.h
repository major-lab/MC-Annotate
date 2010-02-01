/*
 * CycleStatsEntry.h
 *
 *  Created on: Dec 16, 2009
 *      Author: blanchmf
 */

#ifndef CYCLESTATSENTRY_H_
#define CYCLESTATSENTRY_H_

#include <map>
#include <set>
#include <string>

// Prototype
namespace annotate {
	class CycleInfo;
};
class CycleStatsEntry;

class lessCycleStatsEntry
{
public:
	bool operator()(const CycleStatsEntry& aLeft, const CycleStatsEntry& aRight) const;
};

/**
 * An instance of statistics on a given cycle.
 */
class CycleStatsEntry
{
public:
	// LIFECYLE ----------------------------------------------------------------
	CycleStatsEntry(
		const std::string& astrProfile,
		const std::string& astrSequence);
	~CycleStatsEntry();

	// ACCESS ------------------------------------------------------------------
	std::string profile() const { return mstrProfile; }
	std::string sequence() const {	return mstrSequence; }
	unsigned int instances() const { return mInstances.size(); }
	unsigned int interactions() const {	return muiInteractions; }

	// METHODS -----------------------------------------------------------------
	/**
	 * Add an instance of
	 */
	void addInstance(const annotate::CycleInfo& aCycle);
	void addInteraction(
		const annotate::CycleInfo& aThisCycle,
		const annotate::CycleInfo& aOtherCycle);
	unsigned int countInteractingInstances() const;
	unsigned int countInteractions() const;
	std::string toString() const;
private:
	typedef std::set<CycleStatsEntry, lessCycleStatsEntry> cyclestatsentry_set;
	std::string mstrProfile;
	std::string mstrSequence;
	unsigned int muiInteractions;	// Number of interactions of the cycles,
									// only one interaction counted by instance
	std::map<annotate::CycleInfo, cyclestatsentry_set > mInstances;		// Instances of this cycle profile
};

#endif /* CYCLESTATSENTRY_H_ */

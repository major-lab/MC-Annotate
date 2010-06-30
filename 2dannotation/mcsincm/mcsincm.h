/*
 * MCScoreInteractingNCM.h
 *
 *  Created on: June 8, 2010
 *      Author: blanchmf
 */

#ifndef _mcsincm_H_
#define _mcsincm_H_

#include <list>
#include <map>
#include <string>
#include <vector>

#include "AlgorithmExtra.h"
#include "CycleInfo.h"

class MCScoreInteractingNCM
{
public:
	MCScoreInteractingNCM(int argc, char * argv []);

	// METHODS ---------------------------------------------
	void version () const;
	void usage () const;
	void help () const;
private:
	typedef std::vector<annotate::CycleInfo> indexed_cycles;
	typedef std::pair<mccore::ResId, mccore::ResId> resid_pair;
	typedef std::pair<annotate::CycleInfo, annotate::CycleInfo> cycle_pair;
	typedef std::pair<std::pair<float, std::string>, resid_pair> interaction_entry; // (score, face, nuc1, nuc2)

	indexed_cycles mCycles1;
	indexed_cycles mCycles2;

	// std::map<cycle_pair, float> mScores;
	std::map<std::pair<std::string, std::string>, float> mScores;

	float mfAppVersion;
	std::string mstrAppName;

	std::string mstrInteractionsFile;
	std::string mstrCyclesFile1;
	std::string mstrCyclesFile2;

	std::list<interaction_entry> mInteractions;

	std::map<mccore::ResId, std::set<std::string> > mResCycleMap1;
	std::map<mccore::ResId, std::set<std::string> > mResCycleMap2;

	void read_options (int argc, char* argv[]);
	std::list<interaction_entry> readInteractionsFile(const std::string astrFileName) const;

	indexed_cycles readIndexedCycleFile(const std::string& astrFileName) const;

	annotate::CycleInfo::residue_strand getResIds(
		const std::string& aResidues) const;
	annotate::CycleInfo getCycleInfo(
		const std::string& astrIdentifier,
		const std::string& astrProfile,
		const std::string& astrResIds,
		const std::string& astrSequence) const;
	std::string flipSequence(
		unsigned int auiStrand1,
		unsigned int auiStrand2,
		const std::string& astrSequence) const;

	void initializeTables();
	void computeScores();

};
#endif // _mcsincm_H_

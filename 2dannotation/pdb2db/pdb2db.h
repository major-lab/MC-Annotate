/*
 * pdb2db.h
 *
 *  Created on: June 30, 2010
 *      Author: blanchmf
 */

#ifndef _pdb2db_H_
#define _pdb2db_H_

#include <string>
#include <vector>

#include "mccore/Molecule.h"

#include "AnnotationChains.h"
#include "AnnotationStemsLoose.h"

class PDB2DotBracketParams
{
public:
	PDB2DotBracketParams()
	{
		mbOneModel = false;
		muiModelNumber = 0;
		mbBinary = false;
		muiCombinedLayers = 2;
		muiSplitLayers = 0;
		muiMaxPerfectSearch = 24;
	}

	// ATTRIBUTES --------------------------------------------------------------
	bool mbOneModel;
	unsigned int muiModelNumber;  // 1 based vector identifier, 0 means all
	unsigned int muiSplitLayers;
	unsigned int muiCombinedLayers;
	bool mbBinary;
	unsigned int muiMaxPerfectSearch; // Maximum number of stems for exhaustive search

	std::vector<std::string> mFiles;
};


class PDB2DotBracket
{
public:
	// LIFECYCLE -------------------------------------------
	PDB2DotBracket(
		const PDB2DotBracketParams& aParams,
		const std::vector<std::string>& aFiles);
	PDB2DotBracket(
		const PDB2DotBracketParams& aParams,
		const std::string& astrFile,
		std::ostream& aFile);
private:
	typedef std::map<mccore::ResId, char> db_notation;
	typedef std::set<annotate::Stem> stem_set;
	float mfAppVersion;
	std::string mstrAppName;

	PDB2DotBracketParams mParams;
	std::vector<std::string> mFiles;

	// void read_options (int argc, char* argv[]);
	mccore::Molecule* loadFile (const string &filename) const;
	mccore::Molecule* loadStream(const std::ostream& aPDBStream) const;

	void processStream(
		const std::string& astrFile,
		const std::ostream& aFile) const;
	void processFile(const std::string& astrFile) const;
	void processModel(annotate::AnnotateModel& aModel) const;

	std::string getFilePrefix(const std::string& aFileName) const;
	bool pseudoKnots(const std::vector<annotate::Stem>& usedStems, const annotate::Stem& otherStem) const;
	std::string toDotBracketLayers(
		const AnnotationStemsLoose& aStems,
		const std::list<mccore::ResId>& aChain) const;
	std::string toDotBracketCombined(
		const AnnotationStemsLoose& aStems,
		const std::list<mccore::ResId>& aChain) const;
	std::pair<stem_set, stem_set> selectStems(
		const std::vector<annotate::Stem>& aStems) const;
	std::pair<stem_set, stem_set> selectStems(const stem_set& aStems) const;
	std::pair<unsigned int, std::list<std::set<unsigned int> > > selectStemsRecursive(
		const std::set<unsigned int>& aToTest,
		const std::vector<std::set<unsigned int> >& aConflicts,
		const std::vector<annotate::Stem>& aStems) const;
	std::list<stem_set> splitLayers(
		const char acChain,
		const AnnotationStemsLoose& aStems) const;
	std::string toDotBracket(
		const std::list<mccore::ResId>& aChain,
		const std::list<stem_set>& aLayers) const;
	std::string toDotBracket(
		const std::list<mccore::ResId>& aChain,
		const stem_set& aLayer) const;
	void applyStems(
		db_notation& aDBNotation,
		const stem_set& aStems,
		unsigned int auiLevel) const;
	void applyPair(
		db_notation& aDBNotation,
		const annotate::BasePair& aPair,
		const char& acOpen,
		const char& acClose) const;
	db_notation createDotBracket(const std::list<mccore::ResId>& aChain) const;
	std::string getSequence(
		const annotate::AnnotateModel& am,
		const std::list<mccore::ResId>& aChain) const;
	stem_set cutStem(
		const annotate::Stem& aStem,
		const stem_set& aStems) const;
	stem_set cutStems(
		const stem_set& aToCutStem,
		const stem_set& aStems) const;

	stem_set keepChain(
		const char acChain,
		const std::vector<annotate::Stem>& aStems) const;
	stem_set keepChain(
		const char acChain,
		const std::set<annotate::Stem>& aStems) const;
	std::vector<std::set<unsigned int> > computeConflicts(
		const std::vector<annotate::Stem>& aStems) const;
	std::pair<stem_set, stem_set> splitConflicts(const stem_set& aStems) const;

	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > splitImbrication(
		const std::set<annotate::Stem>& aStems) const;

	unsigned int findMaxConflictStem(
		const std::set<unsigned int>& aToTest,
		const std::vector<std::set<unsigned int> >& aConflicts,
		const std::vector<annotate::Stem>& aStems) const;
};

#endif // _pdb2db_H_

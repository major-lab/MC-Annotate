/*
 * pdb2db.h
 *
 *  Created on: June 30, 2010
 *      Author: blanchmf
 */

#ifndef _pdb2db_H_
#define _pdb2db_H_

#include <string>

#include "mccore/Molecule.h"

#include "AnnotationChains.h"
#include "AnnotationStemsLoose.h"

class PDB2DotBracket
{
public:
	PDB2DotBracket(int argc, char * argv []);

	// METHODS ---------------------------------------------
	void version () const;
	void usage () const;
	void help () const;
private:
	typedef std::map<mccore::ResId, char> db_notation;
	typedef std::set<annotate::Stem> stem_set;
	float mfAppVersion;
	std::string mstrAppName;
	std::string mstrPDBFile;

	void read_options (int argc, char* argv[]);
	mccore::Molecule* loadFile (const string &filename);

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
	std::pair<unsigned int, std::list<std::set<unsigned int> > > selectStems(
		const std::set<unsigned int>& aToTest,
		const std::vector<std::set<unsigned int> >& aConflicts,
		const std::vector<annotate::Stem>& aStems) const;
	std::list<stem_set> splitLayers(
		const char acChain,
		const AnnotationStemsLoose& aStems) const;
	std::string toDotBracket(
		const std::list<mccore::ResId>& aChain,
		const std::list<stem_set>& aLayers) const;
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
};

#endif // _pdb2db_H_

/*
 * ModelDotBracket.h
 *
 *  Created on: Dec 20, 2010
 *      Author: blanchmf
 */

#ifndef MODELDOTBRACKET_H_
#define MODELDOTBRACKET_H_

#include "AnnotateModel.h"
#include "AnnotationStemsLoose.h"

namespace annotate
{
	class BasePair;
	class Stem;
};

class ModelDotBracketAnnotator
{
public:
	// LIFECYLE -------------------------------------------
	ModelDotBracketAnnotator(
		bool abCompleteGaps,
		unsigned int auiNbCombinedLayers,
		unsigned int auiNbSplitLayers,
		unsigned int auiMaxPerfectSearch);

	// ACCESSORS ------------------------------------------
	void model(annotate::AnnotateModel& aModel) {mpModel = &aModel;}
	const annotate::AnnotateModel& model() const;

	// METHODS --------------------------------------------
	void process();
private:
	typedef std::pair<unsigned int, unsigned int> index_pair;
	typedef std::list<index_pair> index_pair_list;
	typedef std::map<mccore::ResId, char> db_notation;
	typedef std::set<annotate::Stem> stem_set;

	annotate::AnnotateModel* mpModel;
	bool mbCompleteGaps;
	unsigned int muiNbCombinedLayers;
	unsigned int muiNbSplitLayers;
	unsigned int muiMaxPerfectSearch;
	AnnotationStemsLoose mAnnotationStems;

	std::map<char, index_pair_list> identifyGaps() const;
	std::string insertGapsInString(
		const std::string& aString,
		const index_pair_list& aGaps,
		const char aCharacter) const;
	std::string getSequence(const std::list<mccore::ResId>& aChain) const;
	std::string toDotBracketCombined(
		const std::list<mccore::ResId>& aChain) const;
	std::list<std::string> toDotBracketLayers(
		const std::list<mccore::ResId>& aChain) const;
	bool pseudoKnots(
		const std::vector<annotate::Stem>& usedStems,
		const annotate::Stem& otherStem) const;

	std::pair<stem_set, stem_set> selectStems(
		const std::vector<annotate::Stem>& aStems) const;
	std::pair<stem_set, stem_set> selectStems(
		const stem_set& aStems) const;
	std::pair<unsigned int, std::list<std::set<unsigned int> > > selectStemsRecursive(
		const std::set<unsigned int>& aToTest,
		const std::vector<std::set<unsigned int> >& aConflicts,
		const std::vector<annotate::Stem>& aStems) const;

	std::list<stem_set> splitLayers(const char acChain) const;

	std::list<std::string> toDotBracket(
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

	std::vector<std::set<unsigned int> > computeConflicts(
			const std::vector<annotate::Stem>& aStems) const;
	std::pair<stem_set, stem_set> splitConflicts(const stem_set& aStems) const;

	stem_set cutStem(
		const annotate::Stem& aStem,
		const stem_set& aStems) const;
	stem_set cutStems(
		const stem_set& aToCutStem,
		const stem_set& aStems) const;

	ostream& outputDotBracket(
		std::ostream &aOutputStream,
		const db_notation& aDBNotation) const;

	stem_set keepChain(
		const char acChain,
		const std::vector<annotate::Stem>& aStems) const;
	stem_set keepChain(
		const char acChain,
		const std::set<annotate::Stem>& aStems) const;

	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > splitImbrication(
		const std::set<annotate::Stem>& aStems) const;

	unsigned int findMaxConflictStem(
		const std::set<unsigned int>& aToTest,
		const std::vector<std::set<unsigned int> >& aConflicts,
		const std::vector<annotate::Stem>& aStems) const;

	bool areContiguous(
		const annotate::AnnotateModel& aModel,
		const mccore::Residue& aRes1,
		const mccore::Residue& aRes2) const;
};

#endif /* MODELDOTBRACKET_H_ */

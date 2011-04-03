/*
 * ChainDotBracketAnnotator.h
 *
 *  Created on: Jan 5, 2011
 *      Author: blanchmf
 */

#ifndef CHAINDOTBRACKETANNOTATOR_H_
#define CHAINDOTBRACKETANNOTATOR_H_

#include <libmcannotate/AnnotateModel.h>
#include <libmcannotate/Stem.h>

#include <mccore/ResId.h>

#include <list>

class DBNotation
{
public:
	DBNotation(const std::list<mccore::ResId>& aResidues);
	std::string toString() const;
	void applyPair(
		const annotate::BasePair& aPair,
		const char& acOpen,
		const char& acClose);
private:
	std::map<mccore::ResId, unsigned int> mIndexMapping;
	std::vector<char> mRepresentation;
};

class ChainDotBracketAnnotator
{
public:
	ChainDotBracketAnnotator(
		annotate::AnnotateModel& aModel,
		char acChain,
		const std::list<mccore::ResId>& aResidues,
		bool abConsiderGaps);
	~ChainDotBracketAnnotator() {}

	// METHODS --------------------------------------------
	std::string getSequence(char acGap = 'X') const;
	std::string getDotBracketCombined(
		unsigned int auiNbCombinedLayers,
		unsigned int auiMaxPerfectSearch,
		char acGap = '.') const;
	std::list<std::string> getDotBracketLayers(
		unsigned int auiNbSplitLayers,
		unsigned int auiMaxPerfectSearch,
		char acGap = '.') const;
private:
	typedef std::set<annotate::Stem> stem_set;
	char mcChain;
	std::list<mccore::ResId> mResidues;
	annotate::AnnotateModel* mpModel;
	std::list<std::pair<unsigned int, unsigned int> > mGaps;

	bool areContiguous(
		const mccore::Residue& aRes1,
	  	const mccore::Residue& aRes2) const;

	std::string insertGapsInString(
		const std::string& aString,
		const std::list<std::pair<unsigned int, unsigned int> >& aGaps,
		const char aCharacter) const;
	void identifyGaps();
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > selectStems(
		const std::set<annotate::Stem>& aStems,
		unsigned int auiMaxPerfectSearch) const;
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > selectStems(
		const std::vector<annotate::Stem>& aStems,
		unsigned int auiMaxPerfectSearch) const;

	std::set<annotate::Stem> keepChain(
		const std::vector<annotate::Stem>& aStems) const;
	std::set<annotate::Stem> keepChain(
		const std::set<annotate::Stem>& aStems) const;
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > splitConflicts(
		const std::set<annotate::Stem>& aStems) const;
	std::vector<std::set<unsigned int> > computeConflicts(
		const std::vector<annotate::Stem>& aStems) const;
	bool checkConflict(const annotate::Stem& aStem1, const annotate::Stem& aStem2) const;
	stem_set selectNonConflictedStems(const stem_set& aTarget, const stem_set& aSource) const;
	std::set<annotate::Stem> cutStems(
	 	const std::set<annotate::Stem>& aToCutStems,
	 	const std::set<annotate::Stem>& aStems) const;
	std::set<annotate::Stem> cutStem(
	 	const annotate::Stem& aStem,
	 	const std::set<annotate::Stem>& aStems) const;
	void applyStems(
		DBNotation& aDBNotation,
		const std::set<annotate::Stem>& aStems,
		unsigned int auiLevel) const;
	std::pair<unsigned int, std::list<std::set<unsigned int> > > selectStemsRecursive(
		const std::set<unsigned int>& aToTest,
		const std::vector<std::set<unsigned int> >& aConflicts,
		const std::vector<annotate::Stem>& aStems) const;
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > splitImbrication(
		const std::set<annotate::Stem>& aStems) const;
	std::list<std::set<annotate::Stem> > splitLayers(unsigned int auiMaxPerfectSearch) const;
	std::list<std::string> toDotBracket(
		const std::list<std::set<annotate::Stem> >& aLayers,
		unsigned int auiNbSplitLayers,
		char acGap) const;

	/*
	 * @brief Get the set of residue id contained in a set of stems.
	 * @param aStems : Set of stems from which to get the residue id.
	 * @return Set of residue id found.
	 */
	std::set<mccore::ResId> getResIdSet(const stem_set& aStems) const;
};

#endif /* CHAINDOTBRACKETANNOTATOR_H_ */

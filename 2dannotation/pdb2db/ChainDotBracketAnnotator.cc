/*
 * ChainDotBracketAnnotator.cc
 *
 *  Created on: Jan 5, 2011
 *      Author: blanchmf
 */

#include "ChainDotBracketAnnotator.h"

#include "AnnotationStemsLoose.h"
#include <libmcannotate/Stem.h>

#include <mccorex/AlgorithmExtra.h>

#include <mccore/Pdbstream.h>

#include <cassert>
#include <sstream>

static const unsigned int guiMaxGap = 1000;

static std::string gNotationSymbols[] =
{
	"()ab", // 0
	"[]cd", // 1
	"<>ef", // 2
	"{}gh", // 3
	"+-ij", // 4
	"12kl", // 5
	"34mn", // 6
	"56op", // 7
	"78qr", // 8
	"90st", // 9
};

DBNotation::DBNotation(const std::list<mccore::ResId>& aResidues)
{
	std::list<mccore::ResId>::const_iterator it;
	unsigned int i = 0;
	for(it = aResidues.begin(); it != aResidues.end(); ++ it, ++ i)
	{
		mIndexMapping[*it] = i;
	}

	mRepresentation.resize(aResidues.size(), '.');
}

std::string DBNotation::toString() const
{
	std::ostringstream oss;

	for(unsigned int i = 0; i < mRepresentation.size(); ++ i)
	{
		oss << mRepresentation[i];
	}
	return oss.str();
}

void DBNotation::applyPair(
	const annotate::BasePair& aPair,
	const char& acOpen,
	const char& acClose)
{
	std::map<mccore::ResId, unsigned int>::iterator it1 = mIndexMapping.find(aPair.fResId);
	std::map<mccore::ResId, unsigned int>::iterator it2 = mIndexMapping.find(aPair.rResId);

	if(mRepresentation[it1->second] == '.' && mRepresentation[it2->second] == '.')
	{
		if(it1->second < it2->second)
		{
			mRepresentation[it1->second] = acOpen;
			mRepresentation[it2->second] = acClose;
		}else
		{
			mRepresentation[it1->second] = acClose;
			mRepresentation[it2->second] = acOpen;
		}
	}
}

ChainDotBracketAnnotator::ChainDotBracketAnnotator(
	annotate::AnnotateModel& aModel,
	char acChain,
	const std::list<mccore::ResId>& aResidues,
	bool abConsiderGaps)
: mcChain(acChain)
{
	mResidues = aResidues;
	mpModel = &aModel;
	if(abConsiderGaps)
	{
		identifyGaps();
	}
}

std::string ChainDotBracketAnnotator::getSequence(char acGap) const
{
	assert(0 != mpModel);
	std::ostringstream oss;
	std::list<mccore::ResId>::const_iterator it;
	for(it = mResidues.begin(); it != mResidues.end(); ++ it)
	{
		annotate::AnnotateModel::const_iterator itRes = mpModel->find(*it);
		assert(itRes != mpModel->end());
		oss << mccore::Pdbstream::stringifyResidueType (itRes->getType ());
	}
	std::string strSequence = oss.str();
	strSequence = insertGapsInString(strSequence, mGaps, acGap);
	return strSequence;
}

void ChainDotBracketAnnotator::identifyGaps()
{
	assert(0 != mpModel);

	// Identify the gaps in the sequence, and return the missing ranges
	mGaps.clear();

	try
	{
		bool bBreak = false;
		unsigned int uiPosition = 0;
		annotate::AnnotateModel::const_iterator itPrev = mpModel->begin();
		annotate::AnnotateModel::const_iterator it = itPrev;
		++ it;
		for(; it != mpModel->end() && !bBreak; ++ it)
		{
			uiPosition ++;
			if(it->getResId().getChainId() == itPrev->getResId().getChainId() && !areContiguous(*itPrev, *it))
			{
				if(it->getResId() < itPrev->getResId())
				{
					throw mccore::Exception("Unable to follow sequence to identify gaps.");
				}

				// Nucleotides are not adjacent by atoms
				std::pair<unsigned int, unsigned int> gap;
				gap.first = uiPosition;
				// Find the number of required insertions
				mccore::ResId prevTest = itPrev->getResId();
				prevTest = prevTest + 1;
				unsigned int uiInsertionCount = 1;
				while(!prevTest.areContiguous(it->getResId()))
				{
					if(uiInsertionCount < guiMaxGap)
					{
						prevTest = prevTest + 1;
						uiInsertionCount ++;
					}else
					{
						std::ostringstream oss;
						oss << "Gap at " << gap.first << " too long (" << mpModel->name() << ")";
						oss << " (" << itPrev->getResId() << "," << it->getResId() << ").";
						throw mccore::Exception(oss.str());
					}
				}
				gap.second = uiInsertionCount;
				mGaps.push_back(gap);
			}else if(it->getResId().getChainId() != itPrev->getResId().getChainId())
			{
				uiPosition = 1;
			}
			itPrev = it;
		}
	}catch(const mccore::Exception& e)
	{
		std::cerr << e.GetMessage() << std::endl;
		mGaps.clear();
	}
}

bool ChainDotBracketAnnotator::areContiguous(
  	const mccore::Residue& aRes1,
  	const mccore::Residue& aRes2) const
{
	bool bAreContiguous;
	bAreContiguous = (aRes1.getResId().getChainId() == aRes2.getResId().getChainId());
	if(bAreContiguous)
	{
		if(!aRes1.getResId().areContiguous(aRes2.getResId()))
		{
			// Circumvent a weird side effect of the templates
			mccore::Residue* pRes1 = const_cast<mccore::Residue*>(&aRes1);
			mccore::Residue* pRes2 = const_cast<mccore::Residue*>(&aRes2);
			const mccore::Relation *r = 0;
			try
			{
				r = mpModel->getEdge(pRes1, pRes2);
			}catch(mccore::NoSuchElementException& e){}
			if(0 == r || !(r->is (mccore::PropertyType::pAdjacent)))
			{
				bAreContiguous = false;
			}
		}
	}
	return bAreContiguous;
}

std::string ChainDotBracketAnnotator::insertGapsInString(
	const std::string& aString,
	const std::list<std::pair<unsigned int, unsigned int> >& aGaps,
	const char aCharacter) const
{
	std::string strSequence = aString;
	std::list<std::pair<unsigned int, unsigned int> > gaps = aGaps;
	while(0 < gaps.size())
	{
		std::pair<unsigned int, unsigned int> gap = gaps.back();
		strSequence.insert(gap.first, gap.second, aCharacter);
		gaps.pop_back();
	}
	return strSequence;
}

std::string ChainDotBracketAnnotator::getDotBracketCombined(
	unsigned int auiNbCombinedLayers,
	unsigned int auiMaxPerfectSearch,
	char acGap) const
{
	const AnnotationStemsLoose* pAnnotationStems;
	pAnnotationStems = mpModel->getAnnotation<AnnotationStemsLoose>();
	assert(0 != pAnnotationStems);
	std::list<std::set<annotate::Stem> > layers;

	// Select the non pseudo-knotted stems
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > selectedStems;
	selectedStems.second = keepChain(pAnnotationStems->getStems());
	std::set<annotate::Stem> assignedStems;

	// Stems identified
	while(0 < selectedStems.second.size())
	{
		// Split the conflicting stems from the simple ones
		std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > conflictSplit = splitConflicts(selectedStems.second);
		selectedStems.second = conflictSplit.second;

		// Apply selection on conflicting stems
		selectedStems = selectStems(selectedStems.second, auiMaxPerfectSearch);

		// Keep the selected stems and the unconflicted ones for this layer
		selectedStems.first.insert(conflictSplit.first.begin(), conflictSplit.first.end());
		assignedStems.insert(selectedStems.first.begin(), selectedStems.first.end());
		layers.push_back(selectedStems.first);

		// Remove the already assigned residues from the remaining stems
		unsigned int i = 0;
		std::set<annotate::Stem>::const_iterator it = selectedStems.second.begin();
		for(;it != selectedStems.second.end(); ++ it)
		{
			std::set<annotate::Stem> test;
			test.insert(*it);
			test = cutStems(test, assignedStems);
			++ i;
		}
		selectedStems.second = cutStems(selectedStems.second, assignedStems);
	}

	DBNotation dBrackets(mResidues);
	unsigned int uiLayer = 0;
	std::list<std::set<annotate::Stem> >::const_iterator it;
	for(it = layers.begin(); it != layers.end() && uiLayer < auiNbCombinedLayers; ++ it)
	{
		applyStems(dBrackets, *it, uiLayer);
		uiLayer ++;
	}

	// Output dot-brackets
	return insertGapsInString(dBrackets.toString(), mGaps, acGap);
}

std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > ChainDotBracketAnnotator::selectStems(
	const std::set<annotate::Stem>& aStems,
	unsigned int auiMaxPerfectSearch) const
{
	std::vector<annotate::Stem> stems;
	stems.insert(stems.begin(), aStems.begin(), aStems.end());
	return selectStems(stems, auiMaxPerfectSearch);
}

std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > ChainDotBracketAnnotator::selectStems(
	const std::vector<annotate::Stem>& aStems,
	unsigned int auiMaxPerfectSearch) const
{
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > results;

	if(aStems.size() < auiMaxPerfectSearch)
	{
		std::vector<std::set<unsigned int> > conflicts = computeConflicts(aStems);

		// Set of the stem indices that needs testing
		std::set<unsigned int> toTest;
		for(unsigned int i = 0; i < aStems.size(); ++ i)
		{
			toTest.insert(i);
		}

		std::pair<unsigned int, std::list<std::set<unsigned int> > > selection;
		selection = selectStemsRecursive(toTest, conflicts, aStems);

		for(unsigned int i = 0; i < aStems.size(); ++ i)
		{
			if(selection.second.begin()->find(i) != selection.second.begin()->end())
			{
				results.first.insert(aStems[i]);
			}else
			{
				results.second.insert(aStems[i]);

			}
		}
	}else
	{
		std::set<annotate::Stem> stems;
		stems.insert(aStems.begin(), aStems.end());
		results = splitImbrication(stems);
	}
	return results;
}

std::set<annotate::Stem> ChainDotBracketAnnotator::keepChain(
	const std::vector<annotate::Stem>& aStems) const
{
	std::set<annotate::Stem> keptStems;
	std::set<annotate::Stem> usedStems;
	usedStems.insert(aStems.begin(), aStems.end());
	return keepChain(usedStems);
}

std::set<annotate::Stem> ChainDotBracketAnnotator::keepChain(
	const std::set<annotate::Stem>& aStems) const
{
	std::set<annotate::Stem> keptStems;
	std::set<annotate::Stem>::const_iterator it;
	for(it = aStems.begin(); it != aStems.end(); ++ it)
	{
		char cStemChain = it->basePairs().front().fResId.getChainId();
		if(mcChain == cStemChain)
		{
			keptStems.insert(*it);
		}
	}
	return keptStems;
}

std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > ChainDotBracketAnnotator::splitConflicts(
	const std::set<annotate::Stem>& aStems) const
{
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > returnVal;
	std::vector<std::set<unsigned int> > conflicts;
	std::vector<annotate::Stem> stems;
	stems.insert(stems.begin(), aStems.begin(), aStems.end());

	conflicts = computeConflicts(stems);
	for(unsigned int i = 0; i < stems.size(); ++ i)
	{
		if(0 < conflicts[i].size())
		{
			returnVal.second.insert(stems[i]);
		}else
		{
			returnVal.first.insert(stems[i]);
		}
	}
	return returnVal;
}

std::set<annotate::Stem> ChainDotBracketAnnotator::cutStems(
 	const std::set<annotate::Stem>& aToCutStems,
 	const std::set<annotate::Stem>& aStems) const
 {
	std::set<annotate::Stem> stems;
	std::set<annotate::Stem> cut;
	std::set<annotate::Stem>::const_iterator it;
 	for(it = aToCutStems.begin(); it != aToCutStems.end(); ++ it)
 	{
 		annotate::Stem toCutStem = *it;
 		cut = cutStem(toCutStem, aStems);
 		stems.insert(cut.begin(), cut.end());
 	}
 	return stems;
 }

void ChainDotBracketAnnotator::applyStems(
	DBNotation& aDBNotation,
	const std::set<annotate::Stem>& aStems,
	unsigned int auiLevel) const
{
	std::set< annotate::Stem >::const_iterator it;
	for(it = aStems.begin(); it != aStems.end(); ++ it)
	{
		const std::vector< annotate::BasePair>& pairs = it->basePairs();
		char cOpen = gNotationSymbols[auiLevel][0];
		char cClose = gNotationSymbols[auiLevel][1];
		if(it->getOrientation() == annotate::Stem::ePARALLEL)
		{
			cOpen = gNotationSymbols[auiLevel][2];
			cClose = gNotationSymbols[auiLevel][3];
		}

		std::vector< annotate::BasePair>::const_iterator itPair;
		for(itPair = pairs.begin(); itPair != pairs.end(); ++ itPair)
		{
			aDBNotation.applyPair(*itPair, cOpen, cClose);
		}
	}
}

std::vector<std::set<unsigned int> > ChainDotBracketAnnotator::computeConflicts(
	const std::vector<annotate::Stem>& aStems) const
{
	std::vector<std::set<unsigned int> > conflicts;
	conflicts.resize(aStems.size());
	for(unsigned int i = 0; i < aStems.size(); ++ i)
	{
		for(unsigned int j = 0; j < aStems.size(); ++ j)
		{
			if(i != j)
			{
				if(aStems[i].pseudoKnots(aStems[j]) || aStems[i].overlaps(aStems[j]))
				{
					conflicts[i].insert(j);
				}
			}
		}
	}
	return conflicts;
}

std::pair<unsigned int, std::list<std::set<unsigned int> > > ChainDotBracketAnnotator::selectStemsRecursive(
	const std::set<unsigned int>& aToTest,
	const std::vector<std::set<unsigned int> >& aConflicts,
	const std::vector<annotate::Stem>& aStems) const
{
	std::pair<unsigned int, std::list<std::set<unsigned int> > > returnVal;
	if(0 == aToTest.size())
	{
		std::set<unsigned int> emptySet;
		std::list<std::set<unsigned int> > emptyList;
		returnVal.first = 0;
		returnVal.second.push_back(emptySet);
	}else
	{
		unsigned int uiStem = *aToTest.begin();

		if(0 == aConflicts[uiStem].size())
		{
			// No conflicts
			std::set<unsigned int> subTest = aToTest;
			subTest.erase(uiStem);
			returnVal = selectStemsRecursive(subTest, aConflicts, aStems);
			std::list<std::set<unsigned int> >::iterator it;
			for(it = returnVal.second.begin(); it != returnVal.second.end(); ++ it)
			{
				it->insert(uiStem);
			}
		}else
		{
			std::set<unsigned int> subTest1 = aToTest;
			subTest1.erase(uiStem);

			std::set<unsigned int> subTest2;
			subTest2 = annotate::SetDifference<unsigned int>(subTest1, aConflicts[uiStem]);
			std::pair<unsigned int, std::list<std::set<unsigned int> > > subVal1;
			std::pair<unsigned int, std::list<std::set<unsigned int> > > subVal2;
			subVal1 = selectStemsRecursive(subTest2, aConflicts, aStems);
			subVal2 = selectStemsRecursive(subTest1, aConflicts, aStems);

			// Compute the maximum value of the 2 occurrences
			unsigned int uiMaxValue = std::max(subVal1.first + aStems[uiStem].size(), subVal2.first);
			returnVal.first = uiMaxValue;
			if(uiMaxValue == (subVal1.first + aStems[uiStem].size()))
			{
				returnVal.second.insert(returnVal.second.end(), subVal1.second.begin(), subVal1.second.end());
				std::list<std::set<unsigned int> >::iterator it;
				for(it = returnVal.second.begin(); it != returnVal.second.end(); ++ it)
				{
					it->insert(uiStem);
				}
			}
			if(uiMaxValue == subVal2.first)
			{
				returnVal.second.insert(returnVal.second.end(), subVal2.second.begin(), subVal2.second.end());
			}
		}
	}
	return returnVal;
}

std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > ChainDotBracketAnnotator::splitImbrication(
	const std::set<annotate::Stem>& aStems) const
{
	std::set<annotate::Stem> work = aStems;
	std::set<annotate::Stem> toSelect;
	std::set<annotate::Stem> toReserve;

	while(!work.empty())
	{
		if(1 == work.size())
		{
			// Only one left, keep it
			toSelect.insert(*work.begin());
			work.clear();
		}else
		{
			std::vector<annotate::Stem> stemVector;
			stemVector.insert(stemVector.begin(), work.begin(), work.end());
			std::vector<std::set<unsigned int> > conflicts = computeConflicts(stemVector);
			unsigned int uiMaxStem = 0;
			unsigned int uiMaxCount = 0;
			for(unsigned int i = 0; i < stemVector.size(); ++ i)
			{
				unsigned int uiCount = 0;
				std::set<unsigned int>::const_iterator it;
				for(it = conflicts[i].begin(); it != conflicts[i].end(); ++ it)
				{
					uiCount += stemVector[*it].size();
					if(uiMaxCount < uiCount)
					{
						uiMaxCount = uiCount;
						uiMaxStem = i;
					}
				}
			}
			// Remove max stem
			toReserve.insert(stemVector[uiMaxStem]);
			work.erase(stemVector[uiMaxStem]);
			for(unsigned int j = 0; j < conflicts.size(); ++ j)
			{
				conflicts[j].erase(uiMaxStem);
				if(0 == conflicts[j].size())
				{
					work.erase(stemVector[j]);
					toSelect.insert(stemVector[j]);
				}
			}
		}
	}
	return std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> >(toSelect, toReserve);
}
std::set<annotate::Stem> ChainDotBracketAnnotator::cutStem(
 	const annotate::Stem& aStem,
 	const std::set<annotate::Stem>& aStems) const
 {
 	std::set<annotate::Stem> stems;
 	std::set<annotate::Stem>::const_iterator itStem;

 	std::set<annotate::BasePair> stemPairs;
 	stemPairs.insert(aStem.basePairs().begin(), aStem.basePairs().end());

 	std::set<mccore::ResId> removedResId;
 	for(itStem = aStems.begin(); itStem != aStems.end(); ++ itStem)
 	{
 		std::vector<annotate::BasePair>::const_iterator it;
 		for(it = itStem->basePairs().begin(); it != itStem->basePairs().end(); ++ it)
 		{
 			removedResId.insert(it->fResId);
 			removedResId.insert(it->rResId);
 		}
 	}

 	std::set<annotate::BasePair> remainingPairs;
 	std::set<annotate::BasePair>::const_iterator itPair;
 	for(itPair = stemPairs.begin(); itPair != stemPairs.end(); ++ itPair)
 	{
 		if((removedResId.end() == removedResId.find(itPair->fResId)) && (removedResId.end() == removedResId.find(itPair->rResId)))
 		{
 			remainingPairs.insert(*itPair);
 		}
 	}

 	annotate::Stem stem;
 	for(itPair = remainingPairs.begin(); itPair != remainingPairs.end(); ++ itPair)
 	{
 		if(!stem.continues(*itPair))
 		{
 			stems.insert(stem);
 			stem.clear();
 			stem.push_back(*itPair);
 		}else
 		{
 			stem.push_back(*itPair);
 		}
 	}
 	if(0 < stem.size())
 	{
 		stems.insert(stem);
 		stem.clear();
 	}

 	return stems;
}

std::list<std::string> ChainDotBracketAnnotator::getDotBracketLayers(
	unsigned int auiNbSplitLayers,
	unsigned int auiMaxPerfectSearch,
	char acGap) const
{
	std::list<std::set<annotate::Stem> > layers = splitLayers(auiMaxPerfectSearch);
	return toDotBracket(layers, auiNbSplitLayers, acGap);
}

std::list<std::set<annotate::Stem> > ChainDotBracketAnnotator::splitLayers(unsigned int auiMaxPerfectSearch) const
{
	const AnnotationStemsLoose* pAnnotationStems;
	pAnnotationStems = mpModel->getAnnotation<AnnotationStemsLoose>();
	assert(0 != pAnnotationStems);

	std::list<std::set<annotate::Stem> > layers;
	std::set<annotate::Stem> usedStems;
	usedStems = keepChain(pAnnotationStems->getStems());

	// Select the non pseudo-knotted stems
	std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > selectedStems;
	selectedStems.second = usedStems;

	while(0 < selectedStems.second.size())
	{
		selectedStems = selectStems(selectedStems.second, auiMaxPerfectSearch);
		layers.push_back(selectedStems.first);
	}
	return layers;
}

std::list<std::string> ChainDotBracketAnnotator::toDotBracket(
	const std::list<std::set<annotate::Stem> >& aLayers,
	unsigned int auiNbSplitLayers,
	char acGap) const
{
	std::list<std::string> dotBrackets;

	unsigned int uiLayer = 0;
	std::list<std::set<annotate::Stem> >::const_iterator it;
	for(it = aLayers.begin(); it != aLayers.end() && uiLayer < auiNbSplitLayers; ++ it)
	{
		std::list<mccore::ResId>::const_iterator itResId;
		DBNotation dBrackets(mResidues);

		applyStems(dBrackets, *it, uiLayer);

		uiLayer ++;
		std::string strDB = insertGapsInString(dBrackets.toString(), mGaps, acGap);
		dotBrackets.push_back(strDB);
	}
	return dotBrackets;
}

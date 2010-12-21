/*
 * ModelDotBracket.cc
 *
 *  Created on: Dec 20, 2010
 *      Author: blanchmf
 */

#include "ModelDotBracketAnnotator.h"

#include "AnnotationChains.h"
#include "AnnotationInteractions.h"

#include <mccore/Pdbstream.h>

#include <cassert>
#include <sstream>

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

ModelDotBracketAnnotator::ModelDotBracketAnnotator(
	bool abCompleteGaps,
	unsigned int auiNbCombinedLayers,
	unsigned int auiNbSplitLayers,
	unsigned int auiMaxPerfectSearch)
{
	mpModel = 0;
	mbCompleteGaps = abCompleteGaps;
	muiNbCombinedLayers = auiNbCombinedLayers;
	muiNbSplitLayers = auiNbSplitLayers;
	muiMaxPerfectSearch = auiMaxPerfectSearch;
}

const annotate::AnnotateModel& ModelDotBracketAnnotator::model() const
{
	return *mpModel;
}


void ModelDotBracketAnnotator::process()
{
	assert(0 != mpModel);
	std::map<char, index_pair_list> gaps;
	annotate::AnnotationInteractions annInteractions;
	annotate::AnnotationChains annChains;

	mpModel->addAnnotation(annInteractions);
	mpModel->addAnnotation(annChains);
	mpModel->addAnnotation(mAnnotationStems);

	mpModel->keepRNA();
	unsigned char ucRelationMask =
		mccore::Relation::adjacent_mask
		| mccore::Relation::pairing_mask;
	mpModel->annotate(ucRelationMask);

	if(mbCompleteGaps)
	{
		gaps = identifyGaps();
	}

	// Stems to dot bracket
	annotate::AnnotationChains::chain_map::const_iterator itChain;
	for(itChain = annChains.chains().begin(); itChain != annChains.chains().end(); ++ itChain)
	{
		std::cout << '>' << mpModel->name() << ':' << mpModel->id() << ':';
		std::cout << itChain->first;
		std::cout << "|PDBID|MODEL|CHAIN|SEQUENCE" << std::endl;
		std::string strSequence = getSequence(itChain->second);
		if(mbCompleteGaps)
		{
			strSequence = insertGapsInString(strSequence, gaps[itChain->first], 'X');
		}
		std::cout << strSequence << std::endl;
		std::string strDotBrackets;
		if(0 < muiNbCombinedLayers)
		{
			strDotBrackets = toDotBracketCombined(itChain->second);
			if(mbCompleteGaps)
			{
				strDotBrackets = insertGapsInString(strDotBrackets, gaps[itChain->first], '.');
			}
			std::cout << strDotBrackets << std::endl;
		}
		if(0 < muiNbSplitLayers)
		{
			list<std::string> multiDBs;
			multiDBs = toDotBracketLayers(itChain->second);
			list<std::string>::const_iterator itLayer = multiDBs.begin();
			for(; itLayer != multiDBs.end(); ++ itLayer)
			{
				strDotBrackets = *itLayer;
				if(mbCompleteGaps)
				{
					strDotBrackets = insertGapsInString(strDotBrackets, gaps[itChain->first], '.');
				}
				std::cout << strDotBrackets << std::endl;
			}

		}
	}
}

std::map<char, ModelDotBracketAnnotator::index_pair_list> ModelDotBracketAnnotator::identifyGaps() const
{
	assert(0 != mpModel);

	// Identify the gaps in the sequence, and return the missing ranges
	std::map<char, index_pair_list> gaps;

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
			if(it->getResId().getChainId() == itPrev->getResId().getChainId() && !areContiguous(*mpModel, *itPrev, *it))
			{
				if(it->getResId() < itPrev->getResId())
				{
					throw mccore::Exception("Unable to follow sequence to identify gaps.");
				}

				// Nucleotides are not adjacent by atoms
				index_pair gap;
				gap.first = uiPosition;
				// Find the number of required insertions
				mccore::ResId prevTest = itPrev->getResId();
				prevTest = prevTest + 1;
				unsigned int uiInsertionCount = 1;
				while(!prevTest.areContiguous(it->getResId()))
				{
					if(uiInsertionCount < 1000)
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
				std::map<char, index_pair_list>::iterator itGap;
				itGap = gaps.find(it->getResId().getChainId());
				if(itGap != gaps.end())
				{
					itGap->second.push_back(gap);
				}else
				{
					gaps[it->getResId().getChainId()] = index_pair_list();
					gaps[it->getResId().getChainId()].push_back(gap);
				}
			}else if(it->getResId().getChainId() != itPrev->getResId().getChainId())
			{
				uiPosition = 1;
			}
			itPrev = it;
		}
	}catch(const mccore::Exception& e)
	{
		std::cerr << e.GetMessage() << std::endl;
		gaps.clear();
	}
	return gaps;
}

std::string ModelDotBracketAnnotator::insertGapsInString(
	const std::string& aString,
	const index_pair_list& aGaps,
	const char aCharacter) const
{
	std::string strSequence = aString;
	index_pair_list gaps = aGaps;
	while(0 < gaps.size())
	{
		index_pair gap = gaps.back();
		strSequence.insert(gap.first, gap.second, aCharacter);
		gaps.pop_back();
	}
	return strSequence;
}

std::string ModelDotBracketAnnotator::getSequence(
	const std::list<mccore::ResId>& aChain) const
{
	assert(0 != mpModel);
	std::ostringstream oss;
	std::list<mccore::ResId>::const_iterator it = aChain.begin();
	for(; it != aChain.end(); ++ it)
	{
		annotate::AnnotateModel::const_iterator itRes = mpModel->find(*it);
		assert(itRes != mpModel->end());
		oss << mccore::Pdbstream::stringifyResidueType (itRes->getType ());
	}
	return oss.str();
}

std::string ModelDotBracketAnnotator::toDotBracketCombined(
	const std::list<mccore::ResId>& aChain) const
{
	const char cChain = aChain.begin()->getChainId();
	std::ostringstream oss;
	std::list<stem_set> layers;

	// Select the non pseudo-knotted stems
	std::pair<stem_set, stem_set> selectedStems;
	selectedStems.second = keepChain(cChain, mAnnotationStems.getStems());
	stem_set assignedStems;

	// Stems identified
	while(0 < selectedStems.second.size())
	{
		// Split the conflicting stems from the simple ones
		std::pair<stem_set, stem_set> conflictSplit = splitConflicts(selectedStems.second);
		selectedStems.second = conflictSplit.second;

		// Apply selection on conflicting stems
		selectedStems = selectStems(selectedStems.second);

		// Keep the selected stems and the unconflicted ones for this layer
		selectedStems.first.insert(conflictSplit.first.begin(), conflictSplit.first.end());
		assignedStems.insert(selectedStems.first.begin(), selectedStems.first.end());
		layers.push_back(selectedStems.first);

		// Remove the already assigned residues from the remaining stems
		unsigned int i = 0;
		stem_set::const_iterator it = selectedStems.second.begin();
		for(;it != selectedStems.second.end(); ++ it)
		{
			std::set<annotate::Stem> test;
			test.insert(*it);
			test = cutStems(test, assignedStems);
			++ i;
		}
		selectedStems.second = cutStems(selectedStems.second, assignedStems);
	}

	db_notation dBrackets = createDotBracket(aChain);
	unsigned int uiLayer = 0;
	std::list<stem_set>::const_iterator it;
	for(it = layers.begin(); it != layers.end() && uiLayer < muiNbCombinedLayers; ++ it)
	{
		applyStems(dBrackets, *it, uiLayer);
		uiLayer ++;
	}

	// Output dot-brackets
	outputDotBracket(oss, dBrackets);
	return oss.str();
}

std::list<std::string> ModelDotBracketAnnotator::toDotBracketLayers(
	const std::list<mccore::ResId>& aChain) const
{
	const char cChain = aChain.begin()->getChainId();
	std::list<std::set<annotate::Stem> > layers = splitLayers(cChain);
	return toDotBracket(aChain, layers);
}

bool ModelDotBracketAnnotator::pseudoKnots(
	const std::vector<annotate::Stem>& usedStems,
	const annotate::Stem& otherStem) const
{
	bool bPseudoKnots = false;
	std::vector<annotate::Stem>::const_iterator it;
	for(it = usedStems.begin(); it != usedStems.end() && !bPseudoKnots; ++ it)
	{
		bPseudoKnots = it->pseudoKnots(otherStem);
	}
	return bPseudoKnots;
}

std::pair<ModelDotBracketAnnotator::stem_set, ModelDotBracketAnnotator::stem_set> ModelDotBracketAnnotator::selectStems(
	const stem_set& aStems) const
{
	std::vector<annotate::Stem> stems;
	stems.insert(stems.begin(), aStems.begin(), aStems.end());
	return selectStems(stems);
}

std::pair<ModelDotBracketAnnotator::stem_set, ModelDotBracketAnnotator::stem_set> ModelDotBracketAnnotator::selectStems(
	const std::vector<annotate::Stem>& aStems) const
{
	std::pair<stem_set, stem_set> results;

	if(aStems.size() < muiMaxPerfectSearch)
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

std::list<ModelDotBracketAnnotator::stem_set> ModelDotBracketAnnotator::splitLayers(
	const char acChain) const
{
	std::list<stem_set> layers;
	std::set<annotate::Stem> usedStems;
	usedStems = keepChain(acChain, mAnnotationStems.getStems());

	// Select the non pseudo-knotted stems
	std::pair<stem_set, stem_set> selectedStems;
	selectedStems.second = usedStems;

	while(0 < selectedStems.second.size())
	{
		selectedStems = selectStems(selectedStems.second);
		layers.push_back(selectedStems.first);
	}
	return layers;
}

std::list<std::string> ModelDotBracketAnnotator::toDotBracket(
	const std::list<mccore::ResId>& aChain,
	const std::list<std::set<annotate::Stem> >& aLayers) const
{
	std::list<std::string> dotBrackets;

	unsigned int uiLayer = 0;
	std::list<std::set<annotate::Stem> >::const_iterator it;
	for(it = aLayers.begin(); it != aLayers.end() && uiLayer < muiNbSplitLayers; ++ it)
	{
		std::ostringstream oss;
		std::list<mccore::ResId>::const_iterator itResId;
		db_notation dBrackets = createDotBracket(aChain);

		applyStems(dBrackets, *it, uiLayer);

		outputDotBracket(oss, dBrackets);
		uiLayer ++;
		dotBrackets.push_back(oss.str());
	}
	return dotBrackets;
}

std::string ModelDotBracketAnnotator::toDotBracket(
	const std::list<mccore::ResId>& aChain,
	const stem_set& aLayer) const
{
	std::ostringstream oss;

	db_notation dBrackets = createDotBracket(aChain);
	applyStems(dBrackets, aLayer, 0);

	std::map<mccore::ResId, char>::iterator itDB;
	for(itDB = dBrackets.begin(); itDB != dBrackets.end(); ++ itDB)
	{
		oss << itDB->second;
	}
	return oss.str();
}

ModelDotBracketAnnotator::stem_set ModelDotBracketAnnotator::keepChain(
	const char acChain,
	const std::vector<annotate::Stem>& aStems) const
{
	stem_set keptStems;
	std::set<annotate::Stem> usedStems;
	usedStems.insert(aStems.begin(), aStems.end());
	return keepChain(acChain, usedStems);
}

ModelDotBracketAnnotator::stem_set ModelDotBracketAnnotator::keepChain(
	const char acChain,
	const std::set<annotate::Stem>& aStems) const
{
	stem_set keptStems;
	std::set<annotate::Stem>::const_iterator it;
	for(it = aStems.begin(); it != aStems.end(); ++ it)
	{
		char cStemChain = it->basePairs().front().fResId.getChainId();
		if(acChain == cStemChain)
		{
			keptStems.insert(*it);
		}
	}
	return keptStems;
}

void ModelDotBracketAnnotator::applyStems(
	ModelDotBracketAnnotator::db_notation& aDBNotation,
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
			applyPair(aDBNotation, *itPair, cOpen, cClose);
		}
	}
}

void ModelDotBracketAnnotator::applyPair(
	ModelDotBracketAnnotator::db_notation& aDBNotation,
	const annotate::BasePair& aPair,
	const char& acOpen,
	const char& acClose) const
{
	std::map<mccore::ResId, char>::iterator it1 = aDBNotation.find(aPair.fResId);
	std::map<mccore::ResId, char>::iterator it2 = aDBNotation.find(aPair.rResId);

	if(it1->second == '.' && it2->second == '.')
	{
		it1->second = acOpen;
		it2->second = acClose;
	}
}

std::vector<std::set<unsigned int> > ModelDotBracketAnnotator::computeConflicts(
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

std::pair<ModelDotBracketAnnotator::stem_set, ModelDotBracketAnnotator::stem_set> ModelDotBracketAnnotator::splitConflicts(
	const stem_set& aStems) const
{
	std::pair<stem_set, stem_set> returnVal;
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

std::pair<unsigned int, std::list<std::set<unsigned int> > > ModelDotBracketAnnotator::selectStemsRecursive(
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

unsigned int ModelDotBracketAnnotator::findMaxConflictStem(
	const std::set<unsigned int>& aToTest,
	const std::vector<std::set<unsigned int> >& aConflicts,
	const std::vector<annotate::Stem>& aStems) const
{
	unsigned int uiStem = *aToTest.begin();
	unsigned int uiMaxConflicts = aConflicts[uiStem].size();
	std::set<unsigned int>::const_iterator it;
	for(it = aToTest.begin(); it != aToTest.end(); ++ it)
	{
		if(aConflicts[*it].size() > uiMaxConflicts)
		{
			uiMaxConflicts = aConflicts[*it].size();
			uiStem = *it;
		}
	}
	return uiStem;
}

std::pair<std::set<annotate::Stem>, std::set<annotate::Stem> > ModelDotBracketAnnotator::splitImbrication(
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

ModelDotBracketAnnotator::db_notation ModelDotBracketAnnotator::createDotBracket(const std::list<mccore::ResId>& aChain) const
{
	std::map<mccore::ResId, char> dBrackets;

	// Initialize the dot-bracket notation to unstructured
	std::list<mccore::ResId>::const_iterator itResId;
	for(itResId = aChain.begin(); itResId != aChain.end(); ++ itResId)
	{
		dBrackets[*itResId] = '.';
	}
	return dBrackets;
}


 ostream& ModelDotBracketAnnotator::outputDotBracket(
	std::ostream &aOutputStream,
	const db_notation& aDBNotation) const
{
	 // Output dot-brackets
	std::map<mccore::ResId, char>::const_iterator itDB;
	for(itDB = aDBNotation.begin(); itDB != aDBNotation.end(); ++ itDB)
	{
		aOutputStream << itDB->second;
	}
	return aOutputStream;
}

 ModelDotBracketAnnotator::stem_set ModelDotBracketAnnotator::cutStem(
 	const annotate::Stem& aStem,
 	const stem_set& aStems) const
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

 ModelDotBracketAnnotator::stem_set ModelDotBracketAnnotator::cutStems(
 	const stem_set& aToCutStems,
 	const stem_set& aStems) const
 {
 	stem_set stems;
 	stem_set cut;
 	stem_set::const_iterator it;
 	for(it = aToCutStems.begin(); it != aToCutStems.end(); ++ it)
 	{
 		annotate::Stem toCutStem = *it;
 		cut = cutStem(toCutStem, aStems);
 		stems.insert(cut.begin(), cut.end());
 	}
 	return stems;
 }

 bool ModelDotBracketAnnotator::areContiguous(
 	const annotate::AnnotateModel& aModel,
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
 				r = aModel.getEdge(pRes1, pRes2);
 			}catch(mccore::NoSuchElementException& e){}
 			if(0 == r || !(r->is (mccore::PropertyType::pAdjacent)))
 			{
 				bAreContiguous = false;
 			}
 		}
 	}
 	return bAreContiguous;
 }

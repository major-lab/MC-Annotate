/*
 * ModelDotBracket.cc
 *
 *  Created on: Dec 20, 2010
 *      Author: blanchmf
 */

#include "ModelDotBracketAnnotator.h"
#include "ChainDotBracketAnnotator.h"

#include "AnnotationChains.h"
#include "AnnotationInteractions.h"

#include <mccore/Pdbstream.h>

#include <cassert>
#include <sstream>

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

	// Stems to dot bracket
	annotate::AnnotationChains::chain_map::const_iterator itChain;
	for(itChain = annChains.chains().begin(); itChain != annChains.chains().end(); ++ itChain)
	{
		ChainDotBracketAnnotator chainAnnotator(*mpModel, itChain->first, itChain->second, mbCompleteGaps);

		std::cout << '>' << mpModel->name() << ':' << mpModel->id() << ':';
		std::cout << itChain->first;
		std::cout << "|PDBID|MODEL|CHAIN|SEQUENCE" << std::endl;
		std::string strSequence = chainAnnotator.getSequence('X');
		std::cout << strSequence << std::endl;
		std::string strDotBrackets;
		if(0 < muiNbCombinedLayers)
		{
			strDotBrackets = chainAnnotator.getDotBracketCombined(
				muiNbCombinedLayers,
				muiMaxPerfectSearch,
				'.');
			std::cout << strDotBrackets << std::endl;
		}
		if(0 < muiNbSplitLayers)
		{
			list<std::string> multiDBs = chainAnnotator.getDotBracketLayers(
				muiNbSplitLayers,
				muiMaxPerfectSearch,
				'.');
			list<std::string>::const_iterator itLayer = multiDBs.begin();
			for(; itLayer != multiDBs.end(); ++ itLayer)
			{
				strDotBrackets = *itLayer;
				std::cout << strDotBrackets << std::endl;
			}

		}
	}
}

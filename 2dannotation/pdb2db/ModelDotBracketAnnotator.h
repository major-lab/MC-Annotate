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

	annotate::AnnotateModel* mpModel;
	bool mbCompleteGaps;
	unsigned int muiNbCombinedLayers;
	unsigned int muiNbSplitLayers;
	unsigned int muiMaxPerfectSearch;
	AnnotationStemsLoose mAnnotationStems;
};

#endif /* MODELDOTBRACKET_H_ */

/*
 * AnnotationStemsLoose.h
 *
 *  Created on: Jul 5, 2010
 *      Author: blanchmf
 */

#ifndef _annotate_AnnotationStemsLoose_H_
#define _annotate_AnnotationStemsLoose_H_

#include <vector>

#include "Annotation.h"
#include "Stem.h"

class AnnotationStemsLoose : public annotate::Annotation
{
public:
	AnnotationStemsLoose();
	virtual ~AnnotationStemsLoose();

	virtual void update(annotate::AnnotateModel& aModel);
	virtual std::string output() const;

	const std::vector< annotate::Stem >& getStems() const;
	static const std::string& AnnotationName() {return mstrAnnotationName;}
	virtual const std::string& annotationName() {return AnnotationName();}
private:
	static std::string mstrAnnotationName;
	std::vector< annotate::Stem > mStems;
	virtual void clear();

	void getPotentialStems(
		const annotate::AnnotateModel& aModel,
		std::vector<annotate::Stem>& aStems) const;
	std::set<annotate::BasePair> filterOutMultiChainsPairs(
		std::set<annotate::BasePair>& aPairs) const;

	std::string stemName(
		const annotate::AnnotateModel& aModel,
		const annotate::Stem& aStem);

	bool shouldFilterPair(const annotate::BasePair& aPair) const;
	std::set<annotate::BasePair> filterPairs(const std::vector<annotate::BasePair>& aPairs) const;
};

#endif // _annotate_AnnotationStemsLoose_H_

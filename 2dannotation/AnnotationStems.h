#ifndef _annotate_AnnotationStems_H_
#define _annotate_AnnotationStems_H_

#include <vector>

#include "Annotation.h"
#include "Stem.h"

namespace annotate
{
	class AnnotationStems : public Annotation
	{
	public:
		AnnotationStems();
		virtual ~AnnotationStems();

		virtual void update(AnnotateModel& aModel);
		virtual std::string output() const;

		const std::vector< Stem >& getStems() const;
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
	private:
		static std::string mstrAnnotationName;
		std::vector< Stem > mStems;
		virtual void clear();

		bool checkFacesStrict(const mccore::Relation &aRelation) const;

		bool isBestPartner(
			const AnnotateModel& aModel,
			const BasePair& aPair,
			const std::set<BasePair>& aPairs) const;

		void getPotentialStems(
			const AnnotateModel& aModel,
			std::vector<Stem>& aStems) const;
		void removeInvalids(std::vector<Stem>& aStems) const;
		void removeOverlaps(std::vector<Stem>& aStems) const;
		std::set<BasePair> filterOutMultiChainsPairs(
			std::set<BasePair>& aPairs) const;

		std::string stemName(const AnnotateModel& aModel, const Stem& aStem);
	};

}
#endif /*_annotate_AnnotationStems_H_*/

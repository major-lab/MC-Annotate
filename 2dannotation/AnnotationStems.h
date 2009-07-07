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
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		
		const std::vector< Stem >& getStems() const;
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
	private:
		static std::string mstrAnnotationName;
		std::vector< Stem > mStems;
		virtual void clear();
		
		std::set< BasePair > getWWBasePairs(const AnnotateModel& aModel) const;
		bool checkNucleotides(const mccore::Relation &aRelation) const;
		bool checkFacesRelax(const mccore::Relation &aRelation) const;
		bool checkFacesStrict(const mccore::Relation &aRelation) const;
		bool checkOrientation(const mccore::Relation &aRelation) const;
		bool isWatsonCrick(
			const mccore::Relation &aRelation,
			bool bStrict) const;
		bool isBestPartner(
			const AnnotateModel& aModel, 
			const BasePair& aPair, 
			const std::set<BasePair>& aPairs) const;
		
		void getPotentialStems(
			const AnnotateModel& aModel, 
			std::vector<Stem>& aStems) const;
		void removeInvalids(std::vector<Stem>& aStems) const;
		void removeOverlaps(std::vector<Stem>& aStems) const;
	};
	
}
#endif /*_annotate_AnnotationStems_H_*/
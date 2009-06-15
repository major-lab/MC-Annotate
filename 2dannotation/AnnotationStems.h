#ifndef _annotate_AnnotationStems_H_
#define _annotate_AnnotationStems_H_

#include "Annotation.h"
#include "Stem.h"
#include <vector>

namespace annotate
{
	class AnnotationStems : public Annotation
	{
	public:
		AnnotationStems();
		virtual ~AnnotationStems();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		virtual const std::string provides() const;
		
		const std::vector< Stem >& getStems() const;
	private:
		std::vector< Stem > mStems;
		virtual void clear();
		
		std::set< BasePair > getWWBasePairs(const AnnotateModel& aModel) const;
		void getPotentialStems(
			const AnnotateModel& aModel, 
			std::vector<Stem>& aStems) const;
		void removeInvalids(std::vector<Stem>& aStems) const;
		void removeOverlaps(std::vector<Stem>& aStems) const;
	};
	
}
#endif /*_annotate_AnnotationStems_H_*/
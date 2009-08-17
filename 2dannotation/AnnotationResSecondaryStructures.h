#ifndef _annotate_AnnotationResSecondaryStructure_H_
#define _annotate_AnnotationResSecondaryStructure_H_

#include <vector>

#include "Annotation.h"
#include "Loop.h"

namespace annotate
{
	class AnnotationResSecondaryStructures : public Annotation
	{
	public:
		AnnotationResSecondaryStructures();
		virtual ~AnnotationResSecondaryStructures();
		
		virtual void update(AnnotateModel& aModel);		
		virtual std::string output() const;
		
		const std::map<mccore::ResId, const SecondaryStructure*>& getMapping() const 
		{
			return mMapping;
		}
		
		const std::vector< Loop >& getLoops() const;
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
	private:
		static std::string mstrAnnotationName;
		std::map<mccore::ResId, const SecondaryStructure*> mMapping;
		virtual void clear();
	};
}

#endif /*_annotate_AnnotationResSecondaryStructure_H_*/

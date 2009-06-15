#ifndef _annotate_AnnotationCycles_H_
#define _annotate_AnnotationCycles_H_

#include "Annotation.h"
#include <mccore/GraphModel.h>
#include <vector>

namespace annotate
{
	
	typedef mccore::GraphModel Cycle;
	
	class AnnotationCycles : public Annotation
	{
	public:
		AnnotationCycles();
		virtual ~AnnotationCycles();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		virtual const std::string provides() const;
		
		const std::vector< Cycle >& getCycles() const;
	private:
		std::vector< Cycle > mCycles;
		virtual void clear();
	};
	
}
#endif /*_annotate_AnnotationCycles_H_*/

#ifndef _annotate_AnnotationLoops_H_
#define _annotate_AnnotationLoops_H_

#include <vector>

#include "Annotation.h"
#include "Loop.h"

namespace annotate
{
	class Linker;
	class AnnotationLinkers;
	class AnnotationLoops : public Annotation
	{
	public:
		AnnotationLoops();
		virtual ~AnnotationLoops();
		
		// ACCESS ---------------------------------------------------------------
		const std::vector< Loop >& getLoops() const;
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
		std::vector< Loop > getLoops(const std::string& aDescription) const;
		
		// METHODS --------------------------------------------------------------
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		
		// Output methods
		std::string outputLoop(const Loop& aLoop) const;
		
	private:
		static std::string mstrAnnotationName;
		std::vector< Loop > mLoops;
		virtual void clear();
		
		void getIncompleteLoops(
			const AnnotateModel& aModel, 
			std::list<Loop>& aLoops) const;
		bool loopSelfComplete(const Loop& aLoop) const;
			
		
	};
}

#endif /*_annotate_AnnotationLoop_H_*/

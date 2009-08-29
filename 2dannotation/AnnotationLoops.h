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
		virtual void update(AnnotateModel& aModel);
		virtual std::string output() const;

		// Output methods
		std::string outputLoop(const Loop& aLoop) const;

	private:
		typedef std::list<Loop> loop_list;
		static std::string mstrAnnotationName;
		std::vector< Loop > mLoops;
		virtual void clear();

		/*
		 * @brief For each linker identified in the model, create a loop.
		 * @return List of loops corresponding to the linkers.
		 */
		loop_list createLoopsFromLinkers(const AnnotateModel& aModel) const;

		/*
		 * @brief Find loop by looking at their connectivity with one another.
		 * @return List of loops connected to others.
		 */
		loop_list findLoopsByConnectivity(loop_list& aPotentials) const;
	};
}

#endif /*_annotate_AnnotationLoop_H_*/

#ifndef _annotate_AnnotationLoops_H_
#define _annotate_AnnotationLoops_H_

#include <vector>

#include "Annotation.h"
#include "AnnotationLinkers.h"
#include "Loop.h"

namespace annotate
{
	class Linker;
	class AnnotationInteractions;
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

		std::list<Loop> findLoops(
			AnnotateModel& aModel,
			const std::vector<AnnotationLinkers::linker_info>& aLinkers) const;

		std::list<Loop> findLoop(
			const AnnotationInteractions& aInteractions,
			const std::vector<AnnotationLinkers::linker_info>& aLinkers,
			const std::vector<AnnotationLinkers::linker_info>::const_iterator& aIt) const;

		void removeLooses(
			const std::list<Linker>& aUsed,
			std::list<Linker>& aLooses) const;

		std::vector<AnnotationLinkers::linker_info>::const_iterator advance(
				const std::vector<AnnotationLinkers::linker_info>& aLinkers,
				const std::vector<AnnotationLinkers::linker_info>::const_iterator& aIt,
				const LabeledResId& aResId) const;

	};
}

#endif /*_annotate_AnnotationLoop_H_*/

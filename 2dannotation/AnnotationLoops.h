#ifndef _annotate_AnnotationLoops_H_
#define _annotate_AnnotationLoops_H_

#include "Annotation.h"
#include "Loop.h"

#include <vector>

namespace annotate
{
	class Linker;
	class AnnotationLoops : public Annotation
	{
	public:
		AnnotationLoops();
		virtual ~AnnotationLoops();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		virtual const std::string provides() const;
		
		const std::vector< Loop >& getLoops() const;
	private:
		std::vector< Loop > mLoops;
		virtual void clear();
		
		std::map<mccore::ResId, const Linker*> getResidueLinkerMap(
			const AnnotateModel& aModel) const;
    	Linker nextLinker(
			const Linker& aLinker,
			const std::map<mccore::ResId, const Linker*>& aResidueLinkerMap) const;
		void removeLinker(
			std::map<mccore::ResId, const Linker*>& aResidueLinkerMap,
			const Linker& aLinker) const;
		mccore::ResId nextId(
			const Stem& aStem, 
			const StemConnection& aConnection) const;
			
		// Output methods
		void dumpLoop(std::ostringstream& oss, const Loop& aLoop) const;
	};
}

#endif /*_annotate_AnnotationLoop_H_*/

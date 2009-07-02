#ifndef _annotate_AnnotationLoops_H_
#define _annotate_AnnotationLoops_H_

#include "Annotation.h"
#include "Loop.h"

#include <vector>

namespace annotate
{
	class Linker;
	class AnnotationLinkers;
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
		
		void findLoops(std::map<mccore::ResId, const Linker*> residueLinkerMap);
		void findOpenLoops(std::vector< Loop >& aOpenLoops);
		
		std::map<mccore::ResId, const Linker*> getResidueLinkerMap(
			const AnnotationLinkers& aAnnotLinkers) const;
    	Linker nextLinker(
			const Linker& aLinker,
			const std::map<mccore::ResId, const Linker*>& aResidueLinkerMap) const;
		void removeLinker(
			std::map<mccore::ResId, const Linker*>& aResidueLinkerMap,
			const Linker& aLinker) const;
		std::vector<Loop>::iterator getLoopStartingBy(
			const mccore::ResId& aResId, 
			std::vector<Loop>& aLoops) const;
		std::vector<Loop>::iterator getLoopEndingBy(
			const mccore::ResId& aResId, 
			std::vector<Loop>& aLoops) const;
			
			
		// Output methods
		void dumpLoop(std::ostringstream& oss, const Loop& aLoop) const;
		std::string describeLoop(const Loop& aLoop) const;
	};
}

#endif /*_annotate_AnnotationLoop_H_*/

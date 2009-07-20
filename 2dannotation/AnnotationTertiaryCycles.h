#ifndef _annotate_AnnotationTertiaryCycles_H_
#define _annotate_AnnotationTertiaryCycles_H_

#include <vector>

#include "Annotation.h"
#include "BasePair.h"
#include "BaseStack.h"
#include "Cycle.h"

namespace annotate
{	
	class AnnotationTertiaryPairs;
	class AnnotationTertiaryStacks;
	
	class AnnotationTertiaryCycles : public Annotation
	{
	public:
		AnnotationTertiaryCycles(unsigned int auiMaxCycleSize);
		virtual ~AnnotationTertiaryCycles();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		
		const std::set< Cycle >& getCycles() const {return mCycles;}
		std::set< Cycle >& getCycles()  {return mCycles;}
		
		const std::set<Cycle>& getConnections(const Cycle& aCycle) const;
		
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
		
		// TODO : Make this a method of cycle
		void getPairs(const Cycle& aCycle, std::set<BaseInteraction>& aPairs) const;
	private:
		static std::string mstrAnnotationName;
		std::set< Cycle > mCycles;
		std::map< Cycle, std::set <Cycle> > mConnects;
		unsigned int muiMaxCycleSize;
		virtual void clear();
		bool isTertiary(
			const std::set<BaseInteraction>& aCyclePairs, 
			const std::set<BasePair>& a3DPairs) const;
		bool isTertiary(
			const std::set<BaseInteraction>& aCyclePairs, 
			const std::set<BaseStack>& a3DStacks) const;
			
		void updateConnections(const AnnotateModel& aModel);
	};	
}
#endif /*_annotate_AnnotationTertiaryCycles_H_*/

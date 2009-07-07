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
		AnnotationTertiaryCycles();
		virtual ~AnnotationTertiaryCycles();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		
		const std::vector< Cycle >& getCycles() const;
		std::vector< Cycle >& getCycles();
		
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
	private:
		static std::string mstrAnnotationName;
		std::vector< Cycle > mCycles;
		virtual void clear();
		bool isTertiary(
			const std::set<BaseInteraction>& aCyclePairs, 
			const std::set<BasePair>& a3DPairs) const;
		bool isTertiary(
			const std::set<BaseInteraction>& aCyclePairs, 
			const std::set<BaseStack>& a3DStacks) const;
		void getPairs(const Cycle& aCycle, std::set<BaseInteraction>& aPairs) const;
					
		void outputCycle(
			std::ostringstream& oss, 
			const Cycle& aCycle) const;
			
		std::string getCycleDescription(const Cycle& aCycle) const;
	};	
}
#endif /*_annotate_AnnotationTertiaryCycles_H_*/

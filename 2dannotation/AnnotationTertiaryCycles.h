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
		
		const std::list< Cycle >& getCycles() const {return mCycles;}
		std::list< Cycle >& getCycles()  {return mCycles;}
		
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
		
		// TODO : Make this a method of cycle
		void getPairs(const Cycle& aCycle, std::set<BaseInteraction>& aPairs) const;
	private:
		static std::string mstrAnnotationName;
		std::list< Cycle > mCycles;
		virtual void clear();
		bool isTertiary(
			const std::set<BaseInteraction>& aCyclePairs, 
			const std::set<BasePair>& a3DPairs) const;
		bool isTertiary(
			const std::set<BaseInteraction>& aCyclePairs, 
			const std::set<BaseStack>& a3DStacks) const;
			
		std::string getCycleDescription(const Cycle& aCycle) const;
	};	
}
#endif /*_annotate_AnnotationTertiaryCycles_H_*/

#ifndef _annotate_AnnotationTertiaryCycles_H_
#define _annotate_AnnotationTertiaryCycles_H_

#include "Annotation.h"
#include "BasePair.h"
#include "Cycle.h"
#include <vector>
namespace annotate
{	
	class AnnotationTertiaryPairs;
	
	class AnnotationTertiaryCycles : public Annotation
	{
	public:
		AnnotationTertiaryCycles();
		virtual ~AnnotationTertiaryCycles();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		virtual const std::string provides() const;
		
		const std::vector< Cycle >& getCycles() const;
		std::vector< Cycle >& getCycles();
	private:
		typedef std::pair<mccore::ResId, mccore::ResId> resid_pair;
		typedef std::set<resid_pair> resid_pair_set;
		std::vector< Cycle > mCycles;
		virtual void clear();
		bool isTertiary(
			const resid_pair_set& aCyclePairs, 
			const std::vector<BasePair>& a3DPairs) const;
		void getPairs(const Cycle& aCycle, resid_pair_set& aPairs) const;
			
		bool compare_less(const BasePair& aLeft, const resid_pair& aRight) const;
		bool compare_less(const resid_pair& aLeft, const BasePair& aRight) const;
			
		void outputCycle(
			std::ostringstream& oss, 
			const Cycle& aCycle) const;
			
		std::string getCycleDescription(const Cycle& aCycle) const;
	};
	
}
#endif /*_annotate_AnnotationTertiaryCycles_H_*/

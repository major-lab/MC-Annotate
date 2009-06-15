#ifndef _annotate_AnnotationTertiaryCycles_H_
#define _annotate_AnnotationTertiaryCycles_H_

#include "Annotation.h"
#include "BasePair.h"
#include <mccore/GraphModel.h>
#include <vector>
namespace annotate
{
	typedef mccore::GraphModel Cycle;
	
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
	private:
		typedef std::pair<mccore::ResId, mccore::ResId> resid_pair;
		typedef std::vector<resid_pair> resid_pair_vector;
		typedef std::set<resid_pair> resid_pair_set;
		std::vector< Cycle > mCycles;
		virtual void clear();
		bool isTertiary(
			const resid_pair_vector& aCyclePairs, 
			const resid_pair_vector& a3DPairs) const;
		bool isTertiary(
			const resid_pair& aCyclePair,
			const resid_pair_vector& a3DPairs) const;
		resid_pair_vector getPairs(const Cycle& aCycle) const;
		resid_pair_vector getPairs(
			const AnnotationTertiaryPairs& aTPAnnotations) const;
			
		void outputCycle(
			std::ostringstream& oss, 
			const Cycle& aCycle) const;
	};
	
}
#endif /*_annotate_AnnotationTertiaryCycles_H_*/

#ifndef _annotate_AnnotationTertiaryPairs_H_
#define _annotate_AnnotationTertiaryPairs_H_

#include "Annotation.h"
#include "SecondaryStructure.h"
#include "BasePair.h"
#include <mccore/ResId.h>
#include <vector>
#include <map>

namespace annotate
{
	class AnnotationTertiaryPairs : public Annotation
	{
	public:
		AnnotationTertiaryPairs();
		virtual ~AnnotationTertiaryPairs();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		virtual const std::string provides() const;
		
		const std::vector< BasePair >& getPairs() const;
	private:
		typedef std::pair<mccore::ResId, const SecondaryStructure*> struct_association;
		typedef std::multimap<mccore::ResId, const SecondaryStructure*> struct_association_map;
		std::vector< BasePair > mPairs;
		virtual void clear();
		
		void addStemsAssociations(
			const AnnotateModel& aModel, 
			struct_association_map& associationMap);
			
		void addLoopsAssociations(
			const AnnotateModel& aModel, 
			struct_association_map& associationMap);
			
		bool areAdjacent(
			const mccore::ResId& aResId1, 
			const mccore::ResId& aResId2, 
			const struct_association_map& associationMap) const;
	};
	
}

#endif /*_annotate_AnnotationTertiaryPairs_H_*/

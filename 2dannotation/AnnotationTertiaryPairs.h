#ifndef _annotate_AnnotationTertiaryPairs_H_
#define _annotate_AnnotationTertiaryPairs_H_


#include <set>
#include <map>

#include "mccore/ResId.h"

#include "Annotation.h"
#include "SecondaryStructure.h"
#include "BasePair.h"


namespace annotate
{
	class AnnotationTertiaryPairs : public Annotation
	{
	public:
		AnnotationTertiaryPairs();
		virtual ~AnnotationTertiaryPairs();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		
		const std::set< BasePair >& getPairs() const;
		
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
	private:
		static std::string mstrAnnotationName;
		typedef std::pair<mccore::ResId, const SecondaryStructure*> struct_association;
		typedef std::multimap<mccore::ResId, const SecondaryStructure*> struct_association_map;
		std::set< BasePair > mPairs;
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
			
		bool areSameLoop(
			const AnnotateModel& aModel,
			const mccore::ResId& aResId1, 
			const mccore::ResId& aResId2) const;
	};
	
}

#endif /*_annotate_AnnotationTertiaryPairs_H_*/

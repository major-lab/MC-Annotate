#ifndef _annotate_AnnotationTertiaryStacks_H_
#define _annotate_AnnotationTertiaryStacks_H_

#include "Annotation.h"
#include "SecondaryStructure.h"
#include "BaseStack.h"
#include <mccore/ResId.h>
#include <set>
#include <map>

namespace annotate
{
	class AnnotationTertiaryStacks : public Annotation
	{
	public:
		AnnotationTertiaryStacks();
		virtual ~AnnotationTertiaryStacks();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
				
		const std::set< BaseStack >& getStacks() const;
		
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
	private:
		static std::string mstrAnnotationName;
		typedef std::pair<mccore::ResId, const SecondaryStructure*> struct_association;
		typedef std::multimap<mccore::ResId, const SecondaryStructure*> struct_association_map;
		std::set< BaseStack > mStacks;
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

#endif /*_annotate_AnnotationTertiaryStacks_H_*/

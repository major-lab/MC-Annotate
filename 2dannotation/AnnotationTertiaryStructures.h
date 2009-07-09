#ifndef _annotate_AnnotationTertiaryStructures_H_
#define _annotate_AnnotationTertiaryStructures_H_

#include "mccore/ResId.h"

#include "Annotation.h"
#include "BasePair.h"
#include "BaseStack.h"
#include "Cycle.h"

namespace annotate
{
	class AnnotationTertiaryStructures : public Annotation
	{
	public:
		AnnotationTertiaryStructures();
		virtual ~AnnotationTertiaryStructures();
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		
		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}
		
		const std::list<mccore::GraphModel>& getModels() const {return mMolecules;}
		
	private:		
		static std::string mstrAnnotationName;
		
		std::set<BasePair> mSinglePairs;
		std::set<BaseStack> mSingleStacks;
		std::list<std::list<Cycle> > mStructures;
		std::list<mccore::GraphModel> mMolecules;
		const mccore::GraphModel* mpModel;
		
		virtual void clear();
		
		void computeModels(const AnnotateModel& aModel);
		
		std::set<BasePair> difference(
			const std::set<BasePair>& aPairs, 
			const std::set<BaseInteraction>& aCyclePairs) const;
			
		std::set<BaseStack> difference(
			const std::set<BaseStack>& a3DStack, 
			const std::set<BaseInteraction>& aCyclePairs) const;
	};
	
}

#endif /*_annotate_AnnotationTertiaryStructures_H_*/

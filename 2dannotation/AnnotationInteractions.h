#ifndef _annotate_AnnotationInteractions_H_
#define _annotate_AnnotationInteractions_H_

#include "Annotation.h"
#include "BaseInteraction.h"
#include "BasePair.h"
#include "BaseStack.h"
#include "BaseLink.h"
#include <mccore/GraphModel.h>
#include <vector>
#include <list>

namespace annotate
{
	class AnnotateModel;
	/*
	typedef std::pair<mccore::ResId, mccore::ResId> ResIdPair;*/
	
	class LessBaseInteractionPtr
	{
	public:
  		bool operator() (
  			const BaseInteraction* lhs, 
  			const BaseInteraction* rhs) const
  		{
  			return (lhs->fResId < rhs->fResId
			  || (lhs->fResId == rhs->fResId && lhs->rResId < rhs->rResId));
  		}
	};
	
	class AnnotationInteractions : public Annotation
	{
	public:
		AnnotationInteractions();
		virtual ~AnnotationInteractions();
		
		// This annotation doesn't depends on any other
		void update(const mccore::GraphModel& aModel);
		
		virtual void update(const AnnotateModel& aModel);		
		virtual std::string output() const;
		virtual const std::string provides() const;
		
		const std::vector< BasePair >& getPairs() const {return mPairs;}
		const std::vector< BaseStack >& getStacks() const {return mStacks;}
		const std::vector< BaseLink >& getLinks() const {return mLinks;}
		
		// NOTE : Return value should use smart pointers
		std::list<const BaseInteraction*> getInteractions(
			const mccore::ResId ref, 
			const mccore::ResId res) const; 
		
	private:
		typedef std::multiset< BaseInteraction*, LessBaseInteractionPtr> interaction_multiset;
		typedef interaction_multiset::const_iterator interaction_multiset_const_it;
		const mccore::GraphModel* mpModel;
		std::vector< BasePair > mPairs;
    	std::vector< BaseStack > mStacks;
		std::vector< BaseLink > mLinks;
    	std::vector< unsigned int > mMarks;
    	interaction_multiset mInteractions;
		virtual void clear();
		
		void outputPairs(std::ostringstream& oss) const;
		void outputPair (
			std::ostringstream& oss, 
			const BasePair& aBasePair) const;
		void outputStacks(std::ostringstream& oss) const;
	};
}

#endif /*_annotate_AnnotationInteractions_H_*/

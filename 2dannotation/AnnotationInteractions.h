#ifndef _annotate_AnnotationInteractions_H_
#define _annotate_AnnotationInteractions_H_

#include "Annotation.h"
#include "AlgorithmExtra.h"
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

	class AnnotationInteractions : public Annotation
	{
	public:
		// LIFECYCLE -----------------------------------------------------------
		AnnotationInteractions();
		virtual ~AnnotationInteractions();

		// ACCESS --------------------------------------------------------------
		const std::vector< BasePair >& pairs() const {return mPairs;}
		const std::vector< BaseStack >& stacks() const {return mStacks;}
		const std::vector< BaseLink >& links() const {return mLinks;}
		const std::set<BasePair>& pairsWW() const {return mWWPairs;}

		static const std::string& AnnotationName() {return mstrAnnotationName;}
		virtual const std::string& annotationName() {return AnnotationName();}

		// METHODS -------------------------------------------------------------

		// This annotation doesn't depends on any other
		void update(const mccore::GraphModel& aModel);

		virtual void update(AnnotateModel& aModel);
		virtual std::string output() const;

		// NOTE : Return value should use smart pointers
		std::list<const BaseInteraction*> getInteractions(
			const mccore::ResId ref,
			const mccore::ResId res) const;

		// NOTE : Return value should use smart pointers
		std::list<const BaseInteraction*> getInteractions(
			const std::set<mccore::ResId>& aResIds) const;

		/**
		 * @brief Checks if two residues have a link relation between them
		 */
		bool areContiguous(
			const mccore::ResId ref,
			const mccore::ResId res) const;

	private:
		static std::string mstrAnnotationName;
		typedef std::multiset< BaseInteraction*, less_ptr<BaseInteraction> > interaction_multiset;
		typedef interaction_multiset::const_iterator interaction_multiset_const_it;
		const mccore::GraphModel* mpModel; // Note : This should be a smart ptr
		std::vector< BasePair > mPairs;
    	std::vector< BaseStack > mStacks;
		std::vector< BaseLink > mLinks;
		std::set<BasePair> mWWPairs;
    	std::vector< unsigned int > mMarks;
    	interaction_multiset mInteractions;
		virtual void clear();

		std::set<BasePair> getWWBasePairs(const mccore::GraphModel& aModel) const;
		bool checkFacesStrict(const mccore::Relation &aRelation) const;
		bool isWatsonCrick(const mccore::Relation &aRelation, bool bStrict) const;
		bool checkNucleotides(const mccore::Relation &aRelation) const;
		bool checkFacesRelax(const mccore::Relation &aRelation) const;
		bool checkOrientation(const mccore::Relation &aRelation) const;

		void outputPairs(std::ostringstream& oss) const;
		void outputPair (
			std::ostringstream& oss,
			const BasePair& aBasePair) const;
		void outputStacks(std::ostringstream& oss) const;
	};
}

#endif /*_annotate_AnnotationInteractions_H_*/

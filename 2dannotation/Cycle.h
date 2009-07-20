#ifndef _annotate_Cycle_H_
#define _annotate_Cycle_H_

#include <mccore/GraphModel.h>

#include "AlgorithmExtra.h"
#include "BaseInteraction.h"
#include "AnnotationInteractions.h"

namespace annotate
{
	class Cycle
	{
	public:		
		// LIFECYCLE ------------------------------------------------------------
		Cycle(const mccore::GraphModel& aModel, unsigned char aucRelationMask);
		virtual ~Cycle();
		
		// ACCESS ---------------------------------------------------------------
		const mccore::GraphModel& getModel() const;
		const std::list<mccore::ResId>& getResidues() const {return mResidues;}
		
		const std::string& name() const {return mName;}
		void name(const std::string& aName) {mName = aName;}
		
		const std::string& modelName() const {return mModelName;}
		void modelName(const std::string& aName) {mModelName = aName;}
		
		const std::vector<unsigned int>& profile() const {return mProfile;}
		
		// OPERATORS ------------------------------------------------------------
		bool operator <(const Cycle& aCycle) const;
		
		// METHODS --------------------------------------------------------------
		bool shareInteractions(const Cycle& aCycle) const;
		bool isSingleChain() const;
		std::string getSequence() const;	
		
	private:
		std::string mName;
		std::string mModelName;
				
		typedef std::set<const BaseInteraction*, less_ptr<const BaseInteraction> > interactions_set;
		typedef interactions_set::const_iterator interactions_set_iterator;
		
		unsigned char mucRelationMask;
		mccore::GraphModel mModel;
		AnnotationInteractions mInteractionsAnnotation;
		interactions_set mInteractions;
		std::vector<unsigned int> mProfile;
		std::list<mccore::ResId> mResidues;
		
		void update();
		void updateProfile();
		void clear();
	};
}

#endif /*_annotate_Linker_H_*/
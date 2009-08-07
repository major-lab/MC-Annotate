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
		Cycle(
			const AnnotateModel& aModel, 
			const std::set<BaseInteraction>& aInteractions,
			unsigned char aucRelationMask);
		
		virtual ~Cycle();
		
		// ACCESS ---------------------------------------------------------------
		const mccore::GraphModel& getModel() const;
		const std::list<mccore::ResId>& getResidues() const {return mResidues;}
		
		const std::string& name() const {return mName;}
		void name(const std::string& aName) {mName = aName;}
		
		const std::string& modelName() const {return mModelName;}
		void modelName(const std::string& aName) {mModelName = aName;}
		
		const std::vector<unsigned int>& profile() const {return mProfile;}
		
		const std::set<BaseInteraction>& getBaseInteractions() const 
		{
			return mInteractions;
		}
		
		// OPERATORS ------------------------------------------------------------
		bool operator <(const Cycle& aCycle) const;
		
		// METHODS --------------------------------------------------------------
		bool shareInteractions(const Cycle& aCycle) const;
		bool isSingleChain() const;
		std::string getSequence() const 
			throw(mccore::NoSuchElementException);
		unsigned int getNbStrands() const {return mProfile.size();}
		
	private:
		Cycle() {}
		std::string mName;
		std::string mModelName;
				
		typedef std::set<BaseInteraction> interactions_set;
		typedef interactions_set::const_iterator interactions_set_iterator;
		
		unsigned char mucRelationMask;
		mccore::GraphModel mModel;
		interactions_set mInteractions;
		std::vector<unsigned int> mProfile;
		std::list<mccore::ResId> mResidues;
		
		void updateProfile();
		void clear();
		void setInteractions(
			const std::list<const BaseInteraction*>& aInteractions);
		mccore::ResId getMinResId(const interactions_set& aInteractions) const
			throw(mccore::NoSuchElementException);
		std::list<mccore::ResId> getOrderedResidues(
			const interactions_set& aInteractions) const;
	};
}

#endif /*_annotate_Linker_H_*/
#ifndef _annotate_Cycle_H_
#define _annotate_Cycle_H_

#include "AlgorithmExtra.h"
#include "BaseInteraction.h"
#include "AnnotationInteractions.h"

#include <mccore/GraphModel.h>

#include <memory>

namespace annotate
{
	class Cycle
	{
	public:
		typedef std::set<const BaseInteraction*, less_ptr<BaseInteraction> > interactions_set;
		typedef interactions_set::const_iterator interactions_set_iterator;
		
		enum enType
		{
			eLOOSE,					// Not really a cycle as it is not closed
			eLOOP,					// Simple loop, 1 strand with connected ends 
			e2STRANDS_ANTIPARALLEL,	// 2 Strands cycle with antiparallel strands
			e2STRANDS_PARALLEL,		// 2 Strands cycle with parallel strands
			eMULTIBRANCH,			// Multiple branch cycle
		};
				
		// LIFECYCLE ------------------------------------------------------------
		Cycle(const mccore::GraphModel& aModel, unsigned char aucRelationMask);
		Cycle(
			const mccore::GraphModel& aModel, 
			const interactions_set& aInteractions,
			unsigned char aucRelationMask);
		
		virtual ~Cycle();
		
		// ACCESS ---------------------------------------------------------------
		const mccore::GraphModel& getModel() const;
		const std::list<mccore::ResId>& residues() const {return mResidues;}
		
		const std::string& name() const {return mName;}
		void name(const std::string& aName) {mName = aName;}
		
		const std::string& modelName() const {return mModelName;}
		void modelName(const std::string& aName) {mModelName = aName;}
		
		const std::vector<unsigned int>& profile() const {return mProfile;}
		
		const interactions_set& interactions() const {return mInteractions;}
		
		// OPERATORS ------------------------------------------------------------
		bool operator <(const Cycle& aCycle) const;
		
		// METHODS --------------------------------------------------------------
		bool shareInteractions(const Cycle& aCycle) const;
		bool isSingleChain() const;
		
		/**
		 * @brief Checks if the cycle is complete
		 */
		bool isClosed() const;
		std::string getSequence() const 
			throw(mccore::NoSuchElementException);
		bool isParallel() const;
		unsigned int getNbStrands() const {return mProfile.size();}
		enType getType() const;
		std::vector<std::vector<mccore::ResId> > getStrands() const;
		
		std::set<BaseInteraction> getBaseInteractions() const;
		
	private:
		Cycle() {}
		std::string mName;
		std::string mModelName;
		
		unsigned char mucRelationMask;
		mccore::GraphModel mModel;
		interactions_set mInteractions;
		std::vector<unsigned int> mProfile;
		std::list<mccore::ResId> mResidues;
		
		std::set<BasePair> mPairs;
		std::set<BaseLink> mLinks;
		std::set<BaseStack> mStacks;
		
		void updateProfile();
		void clear();
		void setInteractions(
			const std::list<const BaseInteraction*>& aInteractions);
		std::list<mccore::ResId> getOrderedResidues(
			const interactions_set& aInteractions) const;
			
		void clearInteractions();
			
		void orderResidues(
			const interactions_set& aInteractions,
			std::list<mccore::ResId>& aResidues) const;
		bool isClosed(
			const interactions_set& aInteractions,
			std::list<mccore::ResId>& aResidues) const;
			
		bool areResiduesLinked(
			const mccore::ResId& aRes1, 
			const mccore::ResId& aRes2) const;
		bool areResiduesStacked(
			const mccore::ResId& aRes1, 
			const mccore::ResId& aRes2) const;
		bool areResiduesPaired(
			const mccore::ResId& aRes1, 
			const mccore::ResId& aRes2) const;
	};
}

#endif /*_annotate_Linker_H_*/
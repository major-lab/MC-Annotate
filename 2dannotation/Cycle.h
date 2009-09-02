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
		typedef std::set<BaseInteraction> interactions_set;
		typedef interactions_set::const_iterator interactions_set_iterator;

		enum enType
		{
			eLOOSE,					// Not really a cycle as it is not closed
			eLOOP,					// Simple loop, 1 strand with connected ends
			e2STRANDS_ANTIPARALLEL,	// 2 Strands cycle with antiparallel strands
			e2STRANDS_PARALLEL,		// 2 Strands cycle with parallel strands
			e2STRANDS_TRIANGLE,		// Only one nucleotide on one side
			eMULTIBRANCH,			// Multiple branch cycle
		};

		// LIFECYCLE ------------------------------------------------------------
		Cycle(const interactions_set& aInteractions);

		virtual ~Cycle();

		// ACCESS ---------------------------------------------------------------
		const std::list<mccore::ResId>& resIds() const {return mResidues;}

		const std::string& name() const {return mName;}
		void name(const std::string& aName) {mName = aName;}

		const std::string& modelName() const {return mModelName;}
		void modelName(const std::string& aName) {mModelName = aName;}

		const std::vector<unsigned int>& profile() const {return mProfile;}

		// OPERATORS ------------------------------------------------------------
		bool operator <(const Cycle& aCycle) const;

		// METHODS --------------------------------------------------------------
		/**
		 * @brief Get the interactions from the cycle (with the type member set)
		 */
		const interactions_set getInteractions() const;

		bool shareInteractions(const Cycle& aCycle) const;
		bool isSingleChain() const;

		/**
		 * @brief Checks if the cycle is complete
		 */
		bool isClosed() const;
		bool isParallel() const;
		unsigned int getNbStrands() const {return mProfile.size();}
		enType getType() const;
		std::vector<std::vector<mccore::ResId> > getStrands() const;

		std::set<BaseInteraction> getBaseInteractions() const;

	private:
		Cycle() {}
		std::string mName;
		std::string mModelName;

		interactions_set mInteractions;
		std::vector<unsigned int> mProfile;
		std::list<mccore::ResId> mResidues;

		void updateProfile();
		void clear();
		void setInteractions(const interactions_set& aInteractions);
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

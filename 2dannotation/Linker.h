#ifndef _annotate_Linker_H_
#define _annotate_Linker_H_

#include "SecondaryStructure.h"
#include "Stem.h"
#include "LabeledResId.h"
#include <vector>

namespace annotate
{
	class Linker : public SecondaryStructure
	{
	public:
		Linker() {}

		Linker(
			const std::vector<LabeledResId>& aResidues,
			const SecondaryStructure* apStartStruct,
			const SecondaryStructure* apEndStruct);

		virtual ~Linker();

		// ACCESS ---------------------------------------------------------------
		const std::vector<LabeledResId>& residues() const {return mResidues;}
		const SecondaryStructure* start() const {return mpStartStruct;}
		void start(const SecondaryStructure* apStruct) {mpStartStruct = apStruct;}
		const SecondaryStructure* end() const {return mpEndStruct;}
		void end(const SecondaryStructure* apStruct) {mpEndStruct = apStruct;}

		// OPERATORS ------------------------------------------------------------

		bool operator== (const Linker& other) const;
		bool operator!= (const Linker& other) const;
		bool operator< (const Linker& other) const;

		// METHODS --------------------------------------------------------------

		/**
    	 * @brief Verify the two linkers can be connected.
    	 * @details Two linkers can be connected if they both have an extremity
    	 * connecting to the same pair of a stem.
    	 */
		bool connects(const Linker& aLinker) const;

		/**
		 * @details Checks if this structure and the one passed in parameter are
		 * the same.
		 */
		virtual bool isSame(const SecondaryStructure& aStruct) const;

		/**
		 * @brief Get the set of all residue Ids shared between this secondary
		 * structure and others.
		 */
		virtual std::set<LabeledResId> getSharedResIds() const;

		bool isEmpty() const;
		bool isAdjacent(const SecondaryStructure& aStruct) const;
		bool contains(const LabeledResId& aResId) const;

		void order();
		void reverse();

		/**
		 * @brief getBaseInteractions forming the linker.
		 * @details BaseInteractions are unspecialized interactions between
		 * residues.  This effectively return only one interaction between any
		 * pair of interacting residues.  The interactions are unqualified, ( no
		 * pairing, stacking, adjacency, etc... ).
		 */
		std::list<BaseInteraction> getBaseInteractions() const
			throw(mccore::FatalIntLibException);



	protected:
		std::vector<LabeledResId> mResidues;
		const SecondaryStructure* mpStartStruct;
		const SecondaryStructure* mpEndStruct;

		void clear();

		BaseInteraction getBaseInteractionWithConnection(
			const LabeledResId& aConnection) const
			throw(mccore::FatalIntLibException);

		BaseInteraction getBaseInteractionBetweenConnections() const
			throw(mccore::FatalIntLibException);
	};
}

#endif /*_annotate_Linker_H_*/

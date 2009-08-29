#ifndef _annotate_Loop_H_
#define _annotate_Loop_H_

#include "SecondaryStructure.h"
#include "AlgorithmExtra.h"
#include "Linker.h"
#include "BaseLink.h"

#include <memory>
#include <vector>


namespace annotate
{
	class Loop : public SecondaryStructure
	{
	public:
		typedef std::set<BaseInteraction*, less_ptr<BaseInteraction> > interactions_set;
		Loop();
		Loop(const Linker& aLinker);
		Loop(const std::vector<Linker>& aLinkers);
		virtual ~Loop();

		// ACCESS ---------------------------------------------------------------
		const std::vector<Linker>& linkers() const {return mLinkers;}

		// OPERATORS ------------------------------------------------------------
		bool operator ==(const Loop& other) const;

		// METHODS --------------------------------------------------------------

		/**
    	 * @brief Verify that two loops can be connected.
    	 * @details If a linker at one of the extremity of a loop connects with
    	 * a linker at the other extremity of a loop, than they can be
    	 * connected.
    	 */
		bool connects(const Loop& aLoop) const;

		/**
		 * @details Checks if this structure and the one passed in parameter are
		 * the same.
		 */
		virtual bool isSame(const SecondaryStructure& aStruct) const;

		bool isAdjacent(const SecondaryStructure& aStruct) const;
		bool contains(const mccore::ResId& aResId) const;

		/**
		 * @brief Verify that the loop is closed.
		 * @details The loop is considered closed if both ends connects.
		 */
		bool closed() const;

		/**
		 * @brief Verify that the loop is complete, that it is either a closed
		 * loop or that both end are open.
		 */
		bool complete() const;

		/**
    	 * @brief Insure that the 5' end of the loop is the first position of
    	 * the first linker.
       	 */
		void order();
		void append(const Loop& aLoop);
		std::string describe() const;
		void reverse();

		/**
		 * @brief getBaseInteractions forming the loop.
		 * @details BaseInteractions are unspecialized interactions between
		 * residues.  This effectively return only one interaction between any
		 * pair of interacting residues.  The interactions are unqualified, ( no
		 * pairing, stacking, adjacency, etc... ).
		 */
		std::set<BaseInteraction> getBaseInteractions() const;


		/**
		 * @brief Get a pair of set containing the links and pairs describing
		 * the loop.
		 */
		std::pair<std::set<BaseLink>, std::set<BasePair> > getInteractions() const;

		/**
		 * @brief Get the res ids of the perimeter of the loop.
		 */
		std::set<mccore::ResId> getResIds() const;

		bool checkIntegrity() const;

	protected:
		std::vector<Linker> mLinkers;

		void clear();

		void getLinkerInteractions(
			const Linker& aLinker,
			std::set<BaseLink>& aInteractions) const;

		/**
		 * @brief Append other loops linkers to this loop, in particular orders.
		 */
		void append_back(const Loop& aLoop);
		void append_back_reverse(const Loop& aLoop);
		void append_front(const Loop& aLoop);
		void append_front_reverse(const Loop& aLoop);
	};
}

#endif /*_annotate_Loop_H_*/

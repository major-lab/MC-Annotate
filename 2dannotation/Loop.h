#ifndef _annotate_Loop_H_
#define _annotate_Loop_H_

#include "SecondaryStructure.h"
#include "Linker.h"
#include <vector>

namespace annotate
{	
	class Loop : public SecondaryStructure
	{
	public:		
		Loop();
		Loop(const Linker& aLinker);
		Loop(const std::vector<Linker>& aLinkers);
		virtual ~Loop();
		
		// ACCESS ---------------------------------------------------------------
		const std::vector<Linker>& getLinkers() const;
		
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
		
		bool isAdjacent(const SecondaryStructure& aStruct) const;
		bool contains(const mccore::ResId& aResId) const;
				
		/**
		 * @brief Verify that the loop is closed.
		 * @details The loop is considered closed if both ends connects.
		 */
		bool closed() const;
		
		/**
		 * @brief Verify that the loop is opened.
		 * @details The loop is considered opened if both ends are opened.
		 */
		bool opened() const;
				
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
		 * @brief Get the res ids of the perimeter of the loop.
		 */
		std::set<mccore::ResId> getResIds() const;
		
	protected:
		std::vector<Linker> mLinkers;
		
		void clear();
		
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

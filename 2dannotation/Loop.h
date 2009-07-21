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
		 * @brief Verify that the loop is complete.
		 * @details The loop is considered complete if both ends connects, or 
		 * if both ends are open.
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

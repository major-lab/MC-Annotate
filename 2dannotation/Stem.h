#ifndef _annotate_Stem_H_
#define _annotate_Stem_H_

#include <vector>

#include "BasePair.h"

namespace annotate
{
	class Stem
	{		
	public:
  		enum enOrientation 
  		{ 
  			ePARALLEL, 
  			eANTIPARALLEL, 
  			eUNDEFINED 
  		};
		
		// LIFECYCLE ------------------------------------------------------------
		Stem () { meOrientation = eUNDEFINED; }
	    ~Stem () { }

		// OPERATORS ------------------------------------------------------------

		// ACCESS ---------------------------------------------------------------
		enOrientation getOrientation () const { return meOrientation; }
		// void setId (unsigned int val) { id = val; }
		
		const std::vector< BasePair >& basePairs() const { return mBasePairs; }

		// METHODS --------------------------------------------------------------

	   	/**
    	 * @brief Appends a base pair at the end of the stem.
    	 * @details The base pair is added at the end of the stem if it
    	 * continues the current stem. The pair is not inserted otherwise.
    	 */
    	void push_back(const BasePair& aBasePair);
    	
    	/**
    	 * @brief Empties the stem structure.
    	 */
    	void clear();
    	
    	/**
    	 * @brief Return the number of base pairs in the stem.
    	 * @return number of base pairs
    	 */
    	int size() const {return mBasePairs.size();}
    	
    	/**
    	 * @brief Checks if given pair continues the current stem.
    	 * @return true if the pair continues the stem.
    	 */
    	bool continues(const BasePair& aBasePair) const;
    	
    	/**
    	 * @brief Checks if this stem forms a pseudoknot
    	 * with the given stem.
    	 * @details 
    	 * @return true if the stems forms a pseudoknot, 
    	 * false otherwise.
    	 */
    	bool pseudoKnots( const Stem& aStem ) const;
	private:
		/**
    	 * @brief Order base pair on increasing residue id.
	     * @return the ordered base pair.
    	 */
		BasePair orderPair(const BasePair& aBasePair) const;
		static enOrientation pairOrientation(
			const BasePair& aBasePair1, 
			const BasePair& aBasePair2);
		
		enOrientation meOrientation;
		std::vector< BasePair > mBasePairs;
    	 
	};
}

#endif /*_annotate_Stem_H_*/

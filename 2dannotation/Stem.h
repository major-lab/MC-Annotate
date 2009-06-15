#ifndef _annotate_Stem_H_
#define _annotate_Stem_H_

#include "SecondaryStructure.h"
#include "BasePair.h"

#include <vector>

namespace annotate
{
	class Stem : public SecondaryStructure
	{		
	public:
		enum enConnection
		{
			eFIRST_STRAND_FRONT_PAIR,
			eFIRST_STRAND_BACK_PAIR,
			eSECOND_STRAND_FRONT_PAIR,
			eSECOND_STRAND_BACK_PAIR,
			eUNDEFINED_CONNECTION // Should not be present for a valid connection
		};
		
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
		/**
    	 * @brief Checks if the other stem is equal to this one
    	 * @return true if both Stems have the same base pairs
    	 */
    	 bool operator ==(const Stem& other) const;

		// ACCESS ---------------------------------------------------------------
		enOrientation getOrientation () const { return meOrientation; }
		
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
    	 * @brief Check if current Stem is adjacent to given structure.
    	 * @return true if both structures are known to be adjacent.
    	 */
    	bool isAdjacent(const SecondaryStructure& aStruct) const;
    	
    	/**
    	 * @brief Checks if given residue is involved in 
    	 * this stem's structure.
    	 * @return true if the residue is in one of the 
    	 * pair constituting the stem, false otherwise.
    	 */
    	bool contains(const mccore::Residue& aResidue) const;
    	bool contains(const ResId& aResId) const;
    	
    	/**
    	 * @brief Checks if given stem shares residue with 
    	 * this one.
    	 * @return true if some residues are shared.
    	 */
    	bool overlaps(const Stem& aStem) const;
    	
    	/**
    	 * @brief Checks if this stem forms a pseudoknot
    	 * with the given stem.
    	 * @details 
    	 * @return true if the stems forms a pseudoknot, 
    	 * false otherwise.
    	 */
    	bool pseudoKnots( const Stem& aStem ) const;
    	
    	std::vector<const Stem*> pseudoKnots( 
    		std::vector<const Stem*>& aStems ) const;
    	
    	std::vector<const Stem*> notPseudoKnotted(
    		std::vector<const Stem*>& aStems ) const;
    		
    	ResId getResidue(const enConnection& ePoint) const 
    		throw (mccore::NoSuchElementException);
    		
    	enConnection getConnection(const ResId& aId) const;
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
	
	class StemConnection 
	{
	public:		
		StemConnection() 
		{
			mpStem = NULL; 
			meConnection = Stem::eUNDEFINED_CONNECTION;
		}
		
		StemConnection(const Stem& aStem, const Stem::enConnection& aeConnect)
		{
			mpStem = &aStem;
			meConnection = aeConnect;
		}
		
		bool isValid() const 
		{ 
			return (NULL != mpStem && meConnection != Stem::eUNDEFINED_CONNECTION);
		}
		
		ResId getResidue() const { return mpStem->getResidue(meConnection); }
		
		const Stem& getStem() const {return *mpStem;}
		const Stem::enConnection getConnection() const { return meConnection; }
		
		bool operator== (const StemConnection &other) const 
		{
			bool bEqual = (mpStem == other.mpStem) 
				&& (meConnection == other.meConnection);
			return bEqual;
		}
		
		int getDirection() const;
		mccore::ResId nextId() const;
	private:
		const Stem* mpStem;
		Stem::enConnection meConnection;		
	};
}

#endif /*_annotate_Stem_H_*/

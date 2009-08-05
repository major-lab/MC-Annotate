#ifndef _annotate_Linker_H_
#define _annotate_Linker_H_

#include "SecondaryStructure.h"
#include "Stem.h"
#include "StemConnection.h"
#include <vector>

namespace annotate
{	
	class Linker : public SecondaryStructure
	{
	public:	
		Linker();
		
		Linker(
			const std::vector<mccore::ResId>& aResidues, 
			const StemConnection& aStart,
			const StemConnection& aEnd);
		virtual ~Linker();
		
		// ACCESS ---------------------------------------------------------------
		const std::vector<mccore::ResId >& getResidues() const;
		const StemConnection& getStart() const {return mStart;}
		const StemConnection& getEnd() const {return mEnd;}
		
		// OPERATORS ------------------------------------------------------------
		
		bool operator== (const Linker& other) const;
		bool operator!= (const Linker& other) const;
		bool operator< (const Linker& other) const;
		
		// METHODS --------------------------------------------------------------

		/**
    	 * @brief Verify the two linkers can be connected.
    	 * @details Two linkers can be connected if they both have an extemity 
    	 * connecting to the same pair of a stem.
    	 */
		bool connects(const Linker& aLinker) const;
		
		bool isEmpty() const;
		bool isAdjacent(const SecondaryStructure& aStruct) const;
		bool contains(const mccore::ResId& aResId) const;
		
		void order();
		void reverse();
		
		/**
		 * @brief getBaseInteractions forming the linker.
		 * @details BaseInteractions are unspecialized interactions between 
		 * residues.  This effectively return only one interaction between any 
		 * pair of interacting residues.  The interactions are unqualified, ( no 
		 * pairing, stacking, adjacency, etc... ).
		 */
		std::set<BaseInteraction> getBaseInteractions() const
			throw(mccore::FatalIntLibException);
		
	protected:
		typedef std::pair<mccore::GraphModel::label, mccore::ResId> res_info;
		std::vector<mccore::ResId> mResidues;
		StemConnection mStart;
		StemConnection mEnd;
		
		void clear();
		
		BaseInteraction getBaseInteractionWithConnection(
			const StemConnection& aConnection) const 
			throw(mccore::FatalIntLibException);
	
		BaseInteraction getBaseInteractionBetweenConnections() const 
			throw(mccore::FatalIntLibException);
			
		res_info getResInfo(const StemConnection& aConnection) const 
			throw(mccore::NoSuchElementException);
			
		void appendInteractionWithConnection(
			std::set<BaseInteraction>& aInteractionSet,
			const StemConnection& aConnection) const 
			throw(mccore::FatalIntLibException);
			
		void appendInteractionBetweenConnections(
			std::set<BaseInteraction>& aInteractionSet) const 
			throw(mccore::FatalIntLibException);
	};
}

#endif /*_annotate_Linker_H_*/
